package app

import com.sun.net.httpserver.HttpExchange
import com.sun.net.httpserver.HttpServer
import java.net.InetSocketAddress
import java.nio.charset.StandardCharsets

class WebServer(private val svc: AppService, private val webRoot: java.io.File) {
    private val server: HttpServer = HttpServer.create()

    fun start(port: Int) {
        server.bind(InetSocketAddress("127.0.0.1", port), 0)
        server.createContext("/api/", ::api)
        server.createContext("/", ::staticRoot)
        server.executor = java.util.concurrent.Executors.newFixedThreadPool(8)
        server.start()
    }

    fun stop() = server.stop(0)

    private fun readBody(ex: HttpExchange): String =
        ex.requestBody.readBytes().toString(StandardCharsets.UTF_8)

    private fun sendJson(ex: HttpExchange, status: Int, payload: Any?) {
        val bytes = Json.stringify(payload).toByteArray(StandardCharsets.UTF_8)
        ex.responseHeaders.set("Content-Type", "application/json; charset=utf-8")
        ex.sendResponseHeaders(status, bytes.size.toLong())
        ex.responseBody.use { it.write(bytes) }
    }

    private fun sendError(ex: HttpExchange, status: Int, message: String) =
        sendJson(ex, status, mapOf("error" to message))

    private fun api(ex: HttpExchange) {
        try {
            val path = ex.requestURI.path.removePrefix("/api/").trim('/')
            val method = ex.requestMethod
            when {
                method == "GET" && path == "health" ->
                    sendJson(ex, 200, mapOf("ok" to true, "rules" to ANNOTATION_RULES_VERSION,
                        "transcriptSet" to TRANSCRIPT_SET_VERSION, "builds" to svc.referenceBuilds()))

                method == "GET" && path == "transcripts" ->
                    sendJson(ex, 200, mapOf("transcriptSetVersion" to TRANSCRIPT_SET_VERSION,
                        "transcripts" to svc.transcripts()))

                method == "GET" && path == "batches" ->
                    sendJson(ex, 200, mapOf("batches" to svc.batches()))

                method == "GET" && path == "alignment-window" -> {
                    val q = queryParams(ex)
                    sendJson(ex, 200, svc.alignmentWindow(
                        q["build"]!!.single(), q["chrom"]!!.single(),
                        q["pos"]!!.single().toInt(), q["refLen"]!!.single().toInt()))
                }
                method == "GET" && path.startsWith("groups") -> handleGroupsGet(ex, path)

                method == "POST" && path == "import-batch" -> {
                    val req = Json.parseObject(readBody(ex))
                    val receipt = svc.importBatch(
                        req["batchId"].toString(),
                        req["build"].toString(),
                        req["vcf"].toString(),
                    )
                    sendJson(ex, 200, receipt.toMap())
                }

                method == "POST" && path == "reviews" -> {
                    val req = Json.parseObject(readBody(ex))
                    val version = svc.submitReview(
                        (req["groupId"] as Number).toLong(),
                        (req["expectedVersion"] as Number).toLong(),
                        req["reviewer"].toString(),
                        req["state"].toString(),
                        (req["note"] as? String) ?: "",
                    )
                    sendJson(ex, 200, mapOf("version" to version))
                }

                method == "GET" && path == "sample-vcf" -> {
                    ex.responseHeaders.set("Content-Type", "text/plain; charset=utf-8")
                    val bytes = AppResources.read("/fixtures/sample_batch.vcf").toByteArray(Charsets.UTF_8)
                    ex.sendResponseHeaders(200, bytes.size.toLong())
                    ex.responseBody.use { it.write(bytes) }
                }

                method == "GET" && path == "sample-chain" -> {
                    ex.responseHeaders.set("Content-Type", "text/plain; charset=utf-8")
                    val bytes = AppResources.read("/fixtures/chains/hg19_to_hg38.chain.tsv").toByteArray(Charsets.UTF_8)
                    ex.sendResponseHeaders(200, bytes.size.toLong())
                    ex.responseBody.use { it.write(bytes) }
                }

                method == "POST" && path == "import-chain" -> {
                    val req = Json.parseObject(readBody(ex))
                    val impact = svc.importChain(req["chainTsv"].toString())
                    sendJson(ex, 200, impact.toMap())
                }

                else -> sendError(ex, 404, "no such API route: $method /api/$path")
            }
        } catch (t: Throwable) {
            sendError(ex, 400, t.message ?: t.javaClass.simpleName)
        }
    }

    private fun handleGroupsGet(ex: HttpExchange, path: String) {
        if (path == "groups") {
            val build = queryParams(ex)["build"]?.firstOrNull()
            sendJson(ex, 200, mapOf("groups" to svc.listGroups(build).map {
                mapOf(
                    "id" to it.id, "build" to it.build, "chrom" to it.chrom,
                    "canonPos" to it.canonPos, "canonRef" to it.canonRef, "canonAlt" to it.canonAlt,
                    "refMismatch" to it.refMismatch, "placementCount" to it.placementCount,
                    "reviewState" to it.reviewState, "evidenceCount" to it.evidenceCount,
                )
            }))
            return
        }
        val detailMatch = Regex("groups/([0-9]+)").matchEntire(path)
        if (detailMatch != null) {
            val g = svc.getGroup(detailMatch.groupValues[1].toLong())
                ?: return sendError(ex, 404, "group not found")
            sendJson(ex, 200, g.toMap())
            return
        }
        sendError(ex, 404, "not found")
    }

    private fun queryParams(ex: HttpExchange): Map<String, List<String>> {
        val q = ex.requestURI.rawQuery ?: return emptyMap()
        return q.split('&').mapNotNull { token ->
            val i = token.indexOf('=')
            if (i < 0) null else java.net.URLDecoder.decode(token.substring(0, i), "UTF-8") to
                java.net.URLDecoder.decode(token.substring(i + 1), "UTF-8")
        }.groupBy({ it.first }, { it.second })
    }

    private fun staticRoot(ex: HttpExchange) {
        if (ex.requestMethod != "GET" && ex.requestMethod != "HEAD") {
            ex.sendResponseHeaders(405, -1); ex.close(); return
        }
        val uri = ex.requestURI.path.trim('/')
        val rel = if (uri.isEmpty()) "index.html" else uri
        val file = java.io.File(webRoot, rel).canonicalFile
        if (!file.path.startsWith(webRoot.canonicalPath) || !file.isFile) {
            ex.sendResponseHeaders(404, -1); ex.close(); return
        }
        val type = when (file.extension) {
            "html" -> "text/html; charset=utf-8"
            "js" -> "application/javascript; charset=utf-8"
            "css" -> "text/css; charset=utf-8"
            else -> "application/octet-stream"
        }
        ex.responseHeaders.set("Content-Type", type)
        ex.sendResponseHeaders(200, file.length())
        if (ex.requestMethod == "GET") file.inputStream().use { input -> ex.responseBody.use { input.copyTo(it) } }
        else ex.responseBody.close()
    }

    private fun BatchReceipt.toMap() = mapOf(
        "batchId" to batchId, "build" to build, "receivedAt" to receivedAt,
        "records" to records, "alleles" to alleles, "groups" to groups,
        "refMismatches" to refMismatches, "replay" to replay,
    )

    private fun ImportImpact.toMap() = mapOf(
        "fromBuild" to fromBuild, "toBuild" to toBuild, "segments" to segments,
        "newUniqueGroups" to newUniqueGroups, "flippedGroups" to flippedGroups,
        "manyToOneGroups" to manyToOneGroups, "interruptedGroups" to interruptedGroups,
        "unmappedGroups" to unmappedGroups, "warnings" to warnings,
    )

    private fun GroupDetail.toMap(): Map<String, Any?> = mapOf(
        "id" to id, "build" to build, "chrom" to chrom,
        "canonPos" to canonPos, "canonRef" to canonRef, "canonAlt" to canonAlt,
        "refMismatch" to refMismatch, "expectedRef" to expectedRef,
        "placements" to placements.map {
            mapOf("chrom" to it.chrom, "pos" to it.pos, "ref" to it.ref,
                "alt" to it.alt, "slide" to it.slide)
        },
        "normSteps" to normSteps.map {
            mapOf("kind" to it.kind, "detail" to it.detail, "pos" to it.pos,
                "ref" to it.ref, "alt" to it.alt)
        },
        "preferenceRule" to preferenceRule, "rulesVersion" to rulesVersion,
        "frozen" to frozen, "reviewState" to reviewState, "reviewVersion" to reviewVersion,
        "evidence" to evidence.map {
            mapOf("rawRecordId" to it.rawRecordId, "batchId" to it.batchId, "source" to it.source,
                "lineNumber" to it.lineNumber, "chrom" to it.chrom, "pos" to it.pos,
                "ref" to it.ref, "altField" to it.altField, "altIndex" to it.altIndex,
                "role" to it.role, "rawLine" to it.rawLine, "info" to it.info)
        },
        "consequences" to consequences.map {
            val steps = (Json.parse(it.pathJson) as List<*>).map { row ->
                @Suppress("UNCHECKED_CAST") val m = row as Map<String, Any?>
                mapOf("step" to m["step"], "detail" to m["detail"])
            }
            mapOf("txId" to it.txId, "txBuild" to it.txBuild, "status" to it.status,
                "failureType" to it.failureType, "region" to it.region,
                "transcriptPos" to it.transcriptPos, "refCdna" to it.refCdna,
                "altCdna" to it.altCdna, "hgvsgLike" to it.hgvsgLike, "path" to steps)
        },
        "liftovers" to liftovers.map {
            mapOf("toBuild" to it.toBuild, "status" to it.status, "dstChrom" to it.dstChrom,
                "dstPos" to it.dstPos, "dstRef" to it.dstRef, "dstAlt" to it.dstAlt,
                "detail" to it.detail, "candidateCount" to it.candidateCount)
        },
    )
}
