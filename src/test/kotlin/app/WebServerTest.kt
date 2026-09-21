package app

import org.junit.jupiter.api.Assertions.*
import org.junit.jupiter.api.Test
import java.net.ServerSocket
import java.net.URI
import java.net.http.HttpClient
import java.net.http.HttpRequest
import java.net.http.HttpResponse

class WebServerTest {
    private fun freePort() = ServerSocket(0).use { it.localPort }

    private fun request(port: Int, path: String, method: String = "GET", body: String? = null): HttpResponse<String> {
        val builder = HttpRequest.newBuilder(URI.create("http://127.0.0.1:$port$path"))
        if (body != null) {
            builder.header("Content-Type", "application/json")
            builder.method(method, HttpRequest.BodyPublishers.ofString(body))
        } else builder.method(method, HttpRequest.BodyPublishers.noBody())
        return HttpClient.newHttpClient().send(builder.build(), HttpResponse.BodyHandlers.ofString())
    }

    @Test
    fun `api supports import idempotent replay review conflict and static page`() {
        val dir = FixtureHarness.tempDir()
        val refs = dir.resolve("refs").also { java.nio.file.Files.createDirectories(it) }
        FixtureHarness::class.java.getResourceAsStream("/fixtures/refs/hg19_chr1.fasta").use {
            java.nio.file.Files.copy(it, refs.resolve("hg19_chr1.fasta"))
        }
        FixtureHarness::class.java.getResourceAsStream("/fixtures/refs/hg38_chr1.fasta").use {
            java.nio.file.Files.copy(it, refs.resolve("hg38_chr1.fasta"))
        }
        val web = dir.resolve("web").also { java.nio.file.Files.createDirectories(it) }
        FixtureHarness::class.java.getResourceAsStream("/web/index.html").use {
            java.nio.file.Files.copy(it, web.resolve("index.html"))
        }
        val hg19 = FastaIO.load(refs.resolve("hg19_chr1.fasta").toFile(), "hg19")
        val hg38 = FastaIO.load(refs.resolve("hg38_chr1.fasta").toFile(), "hg38")
        val (setVersion, txs) = TranscriptIo.parseTsv(AppResources.read("/fixtures/transcripts/transcripts.tsv"))
        val db = Database.open(dir.resolve("state.sqlite"))
        val svc = AppService(db, mapOf("hg19" to hg19, "hg38" to hg38), txs.groupBy { it.build }, setVersion)

        val port = freePort()
        val server = WebServer(svc, web.toFile())
        server.start(port)
        try {
            val health = Json.parseObject(request(port, "/api/health").body())
            @Suppress("UNCHECKED_CAST")
            assertTrue((health["builds"] as List<Any?>).contains("hg19"))

            val vcf = AppResources.read("/fixtures/sample_batch.vcf")
            val payload = Json.stringify(mapOf("batchId" to "WEB1", "build" to "hg19", "vcf" to vcf))
            val r1 = Json.parseObject(request(port, "/api/import-batch", "POST", payload).body())
            assertEquals(false, r1["replay"])
            val r2 = Json.parseObject(request(port, "/api/import-batch", "POST", payload).body())
            assertEquals(true, r2["replay"])

            val groups = Json.parseObject(request(port, "/api/groups").body())["groups"] as List<*>
            assertTrue(groups.size >= 5)
            val gid = (groups.first() as Map<*, *>)["id"].toString()
            val detail = request(port, "/api/groups/$gid").body()
            assertTrue(detail.contains("normSteps"))

            // review v0 -> v1 succeeds, stale v0 retry must be a 400
            val rev = Json.stringify(mapOf("groupId" to gid.toLong(), "expectedVersion" to 0L,
                "reviewer" to "alice", "state" to "FROZEN", "note" to "ok"))
            assertEquals(200, request(port, "/api/reviews", "POST", rev).statusCode())
            val staleResp = request(port, "/api/reviews", "POST", rev)
            assertEquals(400, staleResp.statusCode())
            assertTrue(staleResp.body().contains("version 1"))

            val page = request(port, "/").body()
            assertTrue(page.contains("Pair-wise GSB"))
        } finally {
            server.stop()
        }
    }
}
