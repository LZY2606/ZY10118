package app.model

/** Annotation rules version: bumped whenever projection/normalization semantics
 *  change. Every stored projection keeps the version that produced it so that
 *  old results remain reproducible after an upgrade. */
const val ANNOTATION_RULES_VERSION = "1.0.0"

/** A reference build (e.g. b37) with its loaded contig sequences. */
data class Build(
    val name: String,
    val contigs: Map<String, String>,
) {
    fun sequence(contig: String): String? = contigs[contig]
    fun length(contig: String): Int? = contigs[contig]?.length
}

/** Parsed transcript annotation for one build. */
data class Transcript(
    val id: String,
    val build: String,
    val contig: String,
    val strand: String,                 // "+" or "-"
    val exonStarts: List<Int>,          // 1-based inclusive
    val exonEnds: List<Int>,
    val cdsStart: Int,
    val cdsEnd: Int,
    val cdsExonStarts: List<Int>,       // explicit coding exon pieces
    val cdsExonEnds: List<Int>,
) {
    val exons: List<Pair<Int, Int>> get() = exonStarts.zip(exonEnds)
    val codingPieces: List<Pair<Int, Int>> get() = cdsExonStarts.zip(cdsExonEnds)
    fun inExon(g: Int): Boolean = exons.any { it.first <= g && g <= it.second }
    fun inCoding(g: Int): Boolean = codingPieces.any { it.first <= g && g <= it.second }
}

/** One parsed VCF line, attached to its importing batch. */
data class VcfRecord(
    val batchId: String,
    val lineNumber: Int,
    val contig: String,
    val pos: Int,
    val ids: List<String>,
    val ref: String,
    val alts: List<String>,
    val raw: String,
)

/** Status of an observation against the reference build. */
enum class RefStatus { MATCH, MISMATCH }

/** A single evidence record: one (contig, pos, ref, alt) allele occurrence.
 *  Parent alleles produce one observation per ALT plus the parent record. */
data class Observation(
    val batchId: String,
    val lineNumber: Int,
    val contig: String,
    val pos: Int,
    val ref: String,
    val alt: String,
    val build: String,
    val refStatus: RefStatus,
    val parentLineNumber: Int?,         // null for original, set for split children
    val altIndex: Int,                  // 0 for the original composite row, 1..n for children
    val childCount: Int,                // number of children the parent produced
    val sourceId: String,
)

/** Result of minimal representation + left alignment. */
data class Normalized(
    val pos: Int,
    val ref: String,
    val alt: String,
    val steps: List<NormStep>,
    val equivalentAnchors: List<Int>,  // all physical anchor positions
    val uniquePlacement: Boolean,
)

/** A single auditable normalization transformation step. */
data class NormStep(
    val op: String,                     // TRIM_RIGHT | TRIM_LEFT | LEFT_ALIGN
    val pos: Int,
    val ref: String,
    val alt: String,
    val note: String,
)

/** Stable canonical key identifying an equivalence group within a build. */
data class CanonicalKey(val build: String, val contig: String, val pos: Int,
                        val ref: String, val alt: String) {
    fun stableId(): String =
        "%s|%s|%d|%s|%s".format(build, contig, pos, ref, alt)
}

/** Consequence computed for one canonical allele against one transcript. */
data class Projection(
    val transcriptId: String,
    val build: String,
    val rulesVersion: String,
    val consequence: String,
    val failure: String?,               // non-null when consequence is PROJECTION_FAILED
    val cdsStart: Int?,
    val cdsEnd: Int?,
    val proteinStart: Int?,
    val proteinEnd: Int?,
    val refCodon: String?,
    val altCodon: String?,
    val refAa: String?,
    val altAa: String?,
    val path: List<String>,             // human-readable calculation path
    val affectedExons: List<Int>,
)

/** One segment of a liftover chain. */
data class MapSegment(
    val srcContig: String, val srcStart: Int, val srcEnd: Int,
    val dstContig: String, val dstStart: Int, val dstEnd: Int,
    val dstStrand: String,
)

/** Result of lifting one canonical allele. */
data class LiftResult(
    val status: String,                 // MAPPED | GAP | COLLISION | OUT_OF_MAP | STRAND_FLIPPED
    val destContig: String?,
    val destPos: Int?,
    val destRef: String?,
    val destAlt: String?,
    val segmentsUsed: List<Int>,
    val note: String,
)
