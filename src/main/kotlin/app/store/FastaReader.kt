package app.store

import app.model.Build
import java.io.File

/** Minimal FASTA reader: plain or gzip-free fixtures, one or more contigs. */
object FastaReader {
    fun read(file: File, buildName: String): Build {
        val contigs = LinkedHashMap<String, StringBuilder>()
        var current: StringBuilder? = null
        file.useLines { lines ->
            for (raw in lines) {
                val line = raw.trim()
                if (line.isEmpty()) continue
                if (line.startsWith(">")) {
                    val name = line.substring(1).substringBefore(' ').substringBefore('\t')
                    current = StringBuilder()
                    contigs[name] = current!!
                } else {
                    current!!.append(line.uppercase())
                }
            }
        }
        require(contigs.isNotEmpty()) { "FASTA ${file.name} contains no contigs" }
        return Build(buildName, contigs.mapValues { it.value.toString() })
    }
}
