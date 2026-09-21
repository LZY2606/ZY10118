package app

/** Parser for the stripped-down VCF dialect used by this application. */
data class VcfLine(
    val lineNumber: Int,
    val chrom: String,
    val pos: Int,
    val ref: String,
    val alts: List<String>,
    val info: String,
    val rawLine: String,
) {
    val source: String get() {
        val token = info.split(';').firstOrNull { it.startsWith("SRC=") }
        return token?.substringAfter('=')?.takeIf { it.isNotBlank() } ?: "unknown"
    }
}

object VcfParser {
    fun parse(content: String): List<VcfLine> {
        val out = mutableListOf<VcfLine>()
        content.lineSequence().forEachIndexed { idx, raw ->
            val line = raw.trim()
            val lineNumber = idx + 1
            if (line.isEmpty() || line.startsWith("##")) return@forEachIndexed
            if (line.startsWith("#CHROM")) return@forEachIndexed
            val f = line.split('\t')
            require(f.size >= 5) { "line $lineNumber: expected at least 5 tab-separated columns" }
            val pos = f[1].toIntOrNull() ?: error("line $lineNumber: POS not an integer: ${f[1]}")
            val altField = f[4]
            require(altField.isNotEmpty() && altField != ".") { "line $lineNumber: ALT must not be empty/'.' (gVCF not supported)" }
            val alts = altField.split(',').map { it.uppercase() }
            alts.forEach { alt ->
                Variant(f[0], pos, f[3].uppercase(), alt)
            }
            out += VcfLine(lineNumber, f[0], pos, f[3].uppercase(), alts, if (f.size >= 8) f[7] else "", line)
        }
        return out
    }
}
