package app

import java.io.File
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.StandardCopyOption

fun main(args: Array<String>) {
    var port = 8080
    var home = Path.of(".gsb-data")
    var seedDemo = false
    var i = 0
    while (i < args.size) {
        when (args[i]) {
            "--port" -> { port = args[i + 1].toInt(); i += 2 }
            "--data" -> { home = Path.of(args[i + 1]); i += 2 }
            "--seed" -> { seedDemo = true; i += 1 }
            else -> {
                System.err.println("unknown argument: ${args[i]}")
                System.err.println("usage: MainKt [--port 5318] [--data .gsb-data] [--seed]")
                kotlin.system.exitProcess(2)
            }
        }
    }
    Files.createDirectories(home)
    val refDir = home.resolve("refs").also(Files::createDirectories)
    val webDir = home.resolve("web").also(Files::createDirectories)

    // Extract bundled, immutable fixtures (refs + web assets) into the data dir.
    listOf("fixtures/refs/hg19_chr1.fasta", "fixtures/refs/hg38_chr1.fasta").forEach { res ->
        extractResource("/$res", refDir.resolve(res.substringAfterLast('/')).toFile())
    }
    listOf("web/index.html").forEach { res ->
        extractResource("/$res", webDir.resolve(res.substringAfterLast('/')).toFile())
    }

    val hg19 = FastaIO.load(refDir.resolve("hg19_chr1.fasta").toFile(), "hg19")
    val hg38 = FastaIO.load(refDir.resolve("hg38_chr1.fasta").toFile(), "hg38")
    val txContent = readResource("/fixtures/transcripts/transcripts.tsv")
    val (txSetVersion, transcripts) = TranscriptIo.parseTsv(txContent)
    require(txSetVersion == TRANSCRIPT_SET_VERSION) {
        "transcript fixture version mismatch: $txSetVersion != $TRANSCRIPT_SET_VERSION"
    }
    val byBuild = transcripts.groupBy { it.build }

    val db = Database.open(home.resolve("state.sqlite"))
    val svc = AppService(db, mapOf("hg19" to hg19, "hg38" to hg38), byBuild, txSetVersion)

    if (seedDemo) {
        val vcf = readResource("/fixtures/sample_batch.vcf")
        val r = svc.importBatch("DEMO-BATCH-0001", "hg19", vcf)
        val chain = readResource("/fixtures/chains/hg19_to_hg38.chain.tsv")
        val impact = svc.importChain(chain)
        println("seeded batch replay=${r.replay} groups=${r.groups}; " +
            "chain UNIQUE=${impact.newUniqueGroups} FLIPPED=${impact.flippedGroups} " +
            "MANY_TO_ONE=${impact.manyToOneGroups} INTERRUPTED=${impact.interruptedGroups} " +
            "UNMAPPED=${impact.unmappedGroups}")
    }

    val server = WebServer(svc, webDir.toFile())
    server.start(port)
    println("pairwise-gsb listening on http://127.0.0.1:$port")
    println("data directory: ${home.toAbsolutePath()} (annotation rules $ANNOTATION_RULES_VERSION)")
}

private fun extractResource(resource: String, target: File) {
    AppResources::class.java.getResourceAsStream(resource).use { inp ->
        require(inp != null) { "missing bundled resource $resource" }
        Files.copy(inp, target.toPath(), StandardCopyOption.REPLACE_EXISTING)
    }
}

private fun readResource(resource: String): String =
    AppResources::class.java.getResourceAsStream(resource)?.bufferedReader(Charsets.UTF_8)?.use { it.readText() }
        ?: error("missing bundled resource $resource")
