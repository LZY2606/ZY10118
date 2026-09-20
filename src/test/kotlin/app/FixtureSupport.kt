package app

import app.model.Build
import app.model.Transcript
import app.store.FastaReader
import app.transcript.TranscriptLoader
import java.io.File

/** Loads the fixed, dependency-free fixtures shipped under resources. */
object FixtureSupport {
    private fun root(): File {
        val direct = File("src/main/resources/fixtures")
        if (direct.isDirectory) return direct
        val url = javaClass.getResource("/fixtures/ref_b37.fa")
        if (url != null && url.protocol == "file") return File(url.toURI()).parentFile
        error("fixture directory not found")
    }
    val dir = root()
    fun b37(): Build = FastaReader.read(File(dir, "ref_b37.fa"), "b37")
    fun b38(): Build = FastaReader.read(File(dir, "ref_b38.fa"), "b38")
    fun transcriptsB37(): List<Transcript> =
        TranscriptLoader.read(File(dir, "transcripts_b37.tsv"), "b37")
    fun transcriptsB38(): List<Transcript> =
        TranscriptLoader.read(File(dir, "transcripts_b38.tsv"), "b38")
    fun seedVcf(): String = File(dir, "seed_b37.vcf").readText()
    fun chain(): String = File(dir, "map_b37_to_b38.tsv").readText()
}
