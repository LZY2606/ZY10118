package app.vcf

import app.model.VcfRecord

/** Parser for the deliberately small VCF subset used by this application.
 *
 * Coordinate conventions follow VCF 4.2: POS is 1-based, REF and ALT include a
 * shared anchor base, symbolic alleles are rejected, and comma-separated ALTs
 * are preserved verbatim for downstream decomposition.
 */
object VcfParser {
    class VcfFormatException(message: String) : IllegalArgumentException(message)

    fun parse(content: String, batchId: String): List<VcfRecord> {
        val records = ArrayList<VcfRecord>()
        var lineNo = 0
        content.split('\n').forEach { rawLine ->
            val line = rawLine.trimEnd('\r')
            lineNo++
            if (line.isBlank() || line.startsWith("#")) return@forEach
            val c = line.split("\t")
            if (c.size < 5) {
                throw VcfFormatException("line $lineNo: expected at least 5 tab-separated columns")
            }
            val pos = c[1].toIntOrNull()
                ?: throw VcfFormatException("line $lineNo: POS is not an integer: ${c[1]}")
            if (pos < 1) throw VcfFormatException("line $lineNo: POS must be >= 1")
            val ref = c[3].uppercase()
            val alts = c[4].split(",").map { it.uppercase() }
            if (ref.isBlank() || alts.any { it.isBlank() }) {
                throw VcfFormatException("line $lineNo: REF/ALT must not be empty")
            }
            // Structural validation only: alphabetic REF/ALT are accepted so
            // that reference-mismatched records (e.g. wrong REF bases) reach the
            // REF verification stage instead of being rejected at parse time.
            if ((ref + alts.joinToString("")).any { !it.isLetter() }) {
                throw VcfFormatException("line $lineNo: alleles must be nucleotide letters")
            }
            if (alts.any { it.startsWith("<") || it == "*" }) {
                throw VcfFormatException("line $lineNo: symbolic/star alleles are not supported")
            }
            val ids = if (c[2] == ".") emptyList() else c[2].split(";")
            records.add(
                VcfRecord(
                    batchId = batchId,
                    lineNumber = lineNo,
                    contig = c[0],
                    pos = pos,
                    ids = ids,
                    ref = ref,
                    alts = alts,
                    raw = line,
                )
            )
        }
        return records
    }
}
