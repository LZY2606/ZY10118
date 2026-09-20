package app.transcript

import app.model.Transcript
import java.io.File

object TranscriptLoader {
    private fun csv(s: String): List<Int> =
        s.split(",").filter { it.isNotBlank() }.map { it.trim().toInt() }

    fun read(file: File, build: String): List<Transcript> {
        val rows = file.readLines().filter { it.isNotBlank() }
        require(rows.size >= 2) { "transcript TSV ${file.name} has no data rows" }
        return rows.drop(1).map { line ->
            val c = line.split("\t")
            require(c.size >= 10) { "transcript row needs 10 columns: $line" }
            Transcript(
                id = c[0],
                build = build,
                contig = c[1],
                strand = c[2],
                exonStarts = csv(c[3]),
                exonEnds = csv(c[4]),
                cdsStart = c[5].toInt(),
                cdsEnd = c[6].toInt(),
                cdsExonStarts = csv(c[7]),
                cdsExonEnds = csv(c[8]),
            )
        }
    }
}
