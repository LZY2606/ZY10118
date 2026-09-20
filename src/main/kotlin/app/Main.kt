package app

import app.service.AppService
import app.store.Database
import app.web.WebServer
import java.io.File

/** Offline variant normalization review service.
 *
 * Usage:
 *   app.MainKt --port 5318 [--db target/state.db] [--no-seed]
 *
 * On first run the fixed fixtures under src/main/resources/fixtures are loaded
 * automatically so `mvn exec:java` demonstrates every feature end to end.
 */
fun main(args: Array<String>) {
    var port = 5318
    var dbPath = "target/variant-normalizer.db"
    var seed = true
    var reset = false
    var i = 0
    while (i < args.size) {
        when (args[i]) {
            "--port" -> { port = args[++i].toInt() }
            "--db" -> { dbPath = args[++i] }
            "--no-seed" -> seed = false
            "--reset" -> reset = true
            else -> System.err.println("unknown argument: ${args[i]}")
        }
        i++
    }
    if (reset) File(dbPath).delete()

    val db = Database(dbPath)
    val fixtureDir = locateFixtures()
    val svc = AppService(db, fixtureDir)

    if (seed) seedFixtures(svc, fixtureDir)

    val webRoot = locateWebRoot()
    val server = WebServer(svc, webRoot, port)
    server.start()

    Runtime.getRuntime().addShutdownHook(Thread {
        server.stop(); db.close()
    })
    Thread.currentThread().join()
}

private fun seedFixtures(svc: AppService, dir: File) {
    if (svc.queries.listBuilds().any { it["name"] == "b37" }) return
    svc.loadBuild("b37", File(dir, "ref_b37.fa"))
    svc.loadBuild("b38", File(dir, "ref_b38.fa"))
    svc.loadTranscripts("b37", File(dir, "transcripts_b37.tsv"))
    svc.loadTranscripts("b38", File(dir, "transcripts_b38.tsv"))
    val seedVcf = File(dir, "seed_b37.vcf")
    if (seedVcf.isFile) {
        svc.importBatch("SEED-B37-0001", "b37", seedVcf.readText())
    }
    svc.mapImport.importChain("b37", "b38", File(dir, "map_b37_to_b38.tsv").readText())
    println("fixed fixtures seeded (builds b37/b38, batch SEED-B37-0001, b37->b38 chain)")
}

private fun resourceDir(name: String): File? {
    val direct = File("src/main/resources/$name")
    if (direct.isDirectory) return direct
    val url = object {}.javaClass.getResource("/$name/ref_b37.fa")
        ?: object {}.javaClass.getResource("/$name/index.html")
    if (url != null && url.protocol == "file") {
        val f = File(url.toURI()).parentFile
        if (f.isDirectory) return f
    }
    return null
}

private fun locateFixtures(): File =
    resourceDir("fixtures") ?: error("fixture directory not found; run from the project root")

private fun locateWebRoot(): File =
    resourceDir("web") ?: error("web resources not found; run from the project root")
