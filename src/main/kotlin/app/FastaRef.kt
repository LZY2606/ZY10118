package app

import java.io.File

/**
 * Reference genome access. Coordinates are 1-based closed; substring(pos, length)
 * returns `length` bases starting at `pos` and must stay inside the contig.
 */
interface RefGenome {
    val build: String
    fun contigs(): Set<String>
    fun length(chrom: String): Int
    fun substring(chrom: String, pos: Int, length: Int): String
    fun base(chrom: String, pos: Int): Char = substring(chrom, pos, 1).single()
    fun contains(chrom: String, pos: Int, length: Int): Boolean {
        if (chrom !in contigs()) return false
        return pos >= 1 && length >= 1 && pos + length - 1 <= length(chrom)
    }
}

class InMemoryRefGenome(override val build: String, private val seqs: Map<String, String>) : RefGenome {
    constructor(build: String, singleChrom: String) : this(build, mapOf("chr1" to singleChrom))
    override fun contigs() = seqs.keys
    override fun length(chrom: String) = seqs[chrom]?.length ?: error("unknown contig $chrom in build $build")
    override fun substring(chrom: String, pos: Int, length: Int): String {
        val seq = seqs[chrom] ?: error("unknown contig $chrom in build $build")
        require(pos >= 1 && length >= 1 && pos + length - 1 <= seq.length) {
            "out of range: $chrom:$pos len=$length (contig length ${seq.length})"
        }
        return seq.substring(pos - 1, pos - 1 + length).uppercase()
    }
}

object FastaIO {
    fun load(file: File, build: String): InMemoryRefGenome {
        val seqs = LinkedHashMap<String, StringBuilder>()
        var name: String? = null
        file.useLines { lines ->
            for (raw in lines) {
                val line = raw.trim()
                if (line.isEmpty()) continue
                if (line.startsWith(">")) {
                    name = line.removePrefix(">").substringBefore(' ').substringBefore('\t')
                    seqs[name!!] = StringBuilder()
                } else {
                    val sb = seqs[name] ?: error("FASTA sequence before any header in ${file.path}")
                    sb.append(line.uppercase())
                }
            }
        }
        check(seqs.isNotEmpty()) { "no sequences in ${file.path}" }
        return InMemoryRefGenome(build, seqs.mapValues { it.value.toString() })
    }
}
