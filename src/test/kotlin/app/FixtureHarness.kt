package app

import java.io.File
import java.nio.file.Files
import java.nio.file.Path

object FixtureHarness {
    private fun resourceText(name: String): String =
        javaClass.getResourceAsStream(name)!!.bufferedReader(Charsets.UTF_8).use { it.readText() }
    private fun copyResource(name: String, target: Path) {
        javaClass.getResourceAsStream(name).use { inp ->
            requireNotNull(inp) { "missing $name" }
            Files.copy(inp, target)
        }
    }

    fun newService(dir: Path): Triple<AppService, String, String> {
        val refs = dir.resolve("refs").also(Files::createDirectories)
        copyResource("/fixtures/refs/hg19_chr1.fasta", refs.resolve("hg19_chr1.fasta"))
        copyResource("/fixtures/refs/hg38_chr1.fasta", refs.resolve("hg38_chr1.fasta"))
        val hg19 = FastaIO.load(refs.resolve("hg19_chr1.fasta").toFile(), "hg19")
        val hg38 = FastaIO.load(refs.resolve("hg38_chr1.fasta").toFile(), "hg38")
        val (setVersion, txs) = TranscriptIo.parseTsv(resourceText("/fixtures/transcripts/transcripts.tsv"))
        val db = Database.open(dir.resolve("state.sqlite"))
        return Triple(AppService(db, mapOf("hg19" to hg19, "hg38" to hg38), txs.groupBy { it.build }, setVersion),
            resourceText("/fixtures/sample_batch.vcf"),
            resourceText("/fixtures/chains/hg19_to_hg38.chain.tsv"))
    }

    fun tempDir(): Path = Files.createTempDirectory("gsb-test-")
}
