package app.web

import app.service.*
import com.sun.net.httpserver.HttpExchange
import com.sun.net.httpserver.HttpServer
import java.io.File
import java.net.InetSocketAddress
import java.nio.charset.StandardCharsets

class WebServer(
    private val svc: AppService,
    private val webRoot: File,
    port: Int,
) {
    private val server: HttpServer = HttpServer.create(InetSocketAddress("127.0.0.1", port), 0)

    init {
        server.createContext("/") { exchange ->
            try {
                route(exchange)
            } catch (e: Exception) {
                errorResponse(exchange, e)
            }
        }
        server.executor = java.util.concurrent.Executors.newFixedThreadPool(8)
    }

    fun start() {
        server.start()
        println("review UI ready at http://127.0.0.1:${server.address.port}")
    }

    fun stop() = server.stop(0)

    private fun route(ex: HttpExchange) {
        val path = ex.requestURI.path
        val method = ex.requestMethod
        if (method == "GET" && path == "/api/state") return json(ex, 200, mapOf(
            "builds" to svc.queries.listBuilds(),
            "groups" to svc.queries.listGroups(ex.requestURI.query?.let { queryParam(it, "build") }),
            "mismatches" to svc.queries.listMismatches(ex.requestURI.query?.let { queryParam(it, "build") }),
            "rulesVersion" to app.model.ANNOTATION_RULES_VERSION,
        ))
        if (method == "GET" && path.startsWith("/api/groups/")) {
            val id = path.removePrefix("/api/groups/")
            val detail = svc.queries.groupDetail(id)
            val window = try { svc.queries.alignmentWindow(id) } catch (_: Exception) { null }
            return json(ex, 200, detail + ("alignmentWindow" to window))
        }
        if (method == "GET" && path.startsWith("/api/receipts/")) {
            val id = path.removePrefix("/api/receipts/")
            val receipt = svc.queries.batchReceipt(id)
                ?: throw NotFoundException("no receipt for batch $id")
            return json(ex, 200, receipt)
        }
        if (method == "POST" && path == "/api/reviews") {
            val body = bodyJson(ex)
            val cid = body["canonicalId"] as String
            val result = svc.reviews.submit(
                canonicalId = cid,
                reviewer = body["reviewer"] as String,
                action = body["action"] as String,
                expectedVersion = (body["expectedVersion"] as Number).toInt(),
                anchorOverride = (body["anchorOverride"] as Number?)?.toInt(),
                note = body["note"] as String?,
            )
            return json(ex, 200, result)
        }
        if (method == "POST" && path == "/api/batches") {
            val body = bodyJson(ex)
            val receipt = svc.importBatch(
                body["batchId"] as String,
                body["build"] as String,
                body["vcf"] as String,
            )
            return json(ex, 200, receipt.toMap())
        }
        if (method == "POST" && path == "/api/liftover") {
            val body = bodyJson(ex)
            val report = svc.mapImport.importChain(
                body["sourceBuild"] as String,
                body["destBuild"] as String,
                body["chainTsv"] as String,
            )
            return json(ex, 200, mapOf(
                "sourceBuild" to report.sourceBuild,
                "destBuild" to report.destBuild,
                "total" to report.total,
                "summary" to report.summary,
                "items" to report.items.map {
                    mapOf("canonicalId" to it.canonicalId, "status" to it.status,
                        "destContig" to it.destContig, "destPos" to it.destPos,
                        "segments" to it.segments, "note" to it.note)
                },
            ))
        }
        if (method == "GET" && path.startsWith("/api/transcripts/")) {
            val buildName = path.removePrefix("/api/transcripts/")
            return json(ex, 200, mapOf("transcripts" to svc.transcripts(buildName).map {
                mapOf("id" to it.id, "contig" to it.contig, "strand" to it.strand,
                    "exons" to it.exons.map { (a, b) -> mapOf("start" to a, "end" to b) },
                    "codingPieces" to it.codingPieces.map { (a, b) ->
                        mapOf("start" to a, "end" to b) })
            }))
        }
        serveStatic(ex, path)
    }

    private fun queryParam(query: String, key: String): String? =
        query.split("&").map { it.split("=", limit = 2) }
            .firstOrNull { it[0] == key }?.getOrNull(1)

    private fun bodyJson(ex: HttpExchange): Map<String, Any?> {
        val text = ex.requestBody.readBytes().toString(StandardCharsets.UTF_8)
        val parsed = JsonParser.parse(text)
        if (parsed !is Map<*, *>) throw BadRequestException("request body must be a JSON object")
        @Suppress("UNCHECKED_CAST")
        return parsed as Map<String, Any?>
    }

    private fun json(ex: HttpExchange, status: Int, value: Any?) {
        val bytes = Json.stringify(value).toByteArray(StandardCharsets.UTF_8)
        ex.responseHeaders.set("Content-Type", "application/json; charset=utf-8")
        ex.sendResponseHeaders(status, bytes.size.toLong())
        ex.responseBody.use { it.write(bytes) }
    }

    private fun errorResponse(ex: HttpExchange, e: Exception) {
        val (status, code) = when (e) {
            is ConflictException -> 409 to "CONFLICT"
            is BadRequestException -> 400 to "BAD_REQUEST"
            is NotFoundException -> 404 to "NOT_FOUND"
            else -> 500 to "INTERNAL_ERROR"
        }
        if (status == 500) e.printStackTrace()
        json(ex, status, mapOf("error" to code, "message" to (e.message ?: code)))
    }

    private fun serveStatic(ex: HttpExchange, path: String) {
        val rel = if (path == "/") "index.html" else path.trimStart('/')
        val file = File(webRoot, rel).canonicalFile
        if (!file.path.startsWith(webRoot.canonicalPath) || !file.isFile) {
            return json(ex, 404, mapOf("error" to "NOT_FOUND", "message" to rel))
        }
        val type = when (file.extension) {
            "html" -> "text/html; charset=utf-8"
            "js" -> "application/javascript; charset=utf-8"
            "css" -> "text/css; charset=utf-8"
            else -> "application/octet-stream"
        }
        val bytes = file.readBytes()
        ex.responseHeaders.set("Content-Type", type)
        ex.sendResponseHeaders(200, bytes.size.toLong())
        ex.responseBody.use { it.write(bytes) }
    }
}
