package app

const val ANNOTATION_RULES_VERSION = "annot-rules-1.0.0"
const val TRANSCRIPT_SET_VERSION = "tx-fixture-1"

data class Variant(
    val chrom: String,
    val pos: Int,
    val ref: String,
    val alt: String,
) {
    init {
        require(pos >= 1) { "VCF POS must be 1-based, got $pos" }
        require(ref.isNotEmpty() && alt.isNotEmpty()) { "REF/ALT must not be empty" }
        require(ref.all { it in "ACGTNacgtn" } && alt.all { it in "ACGTNacgtn,*" }) {
            "non-DNA character in REF/ALT"
        }
    }
    val key: String get() = "$chrom:$pos:$ref:${alt.uppercase()}"
}

/** One equivalent physical VCF placement of an allele. */
data class Placement(
    val chrom: String,
    val pos: Int,
    val ref: String,
    val alt: String,
    val slide: Int,
) {
    fun variant() = Variant(chrom, pos, ref.uppercase(), alt.uppercase())
}

/** A complete normalization result for one allele within one build. */
data class NormalizedAllele(
    val input: Variant,
    val refMatches: Boolean,
    val expectedRef: String?,
    val steps: List<NormStep>,
    val placements: List<Placement>,
    val canonical: Placement?,
    val preferenceRule: String,
)

data class NormStep(val kind: String, val detail: String, val pos: Int? = null, val ref: String? = null, val alt: String? = null)

data class Transcript(
    val txId: String,
    val gene: String,
    val chrom: String,
    val strand: Char,
    val exons: List<Pair<Int, Int>>,
    val cdsStart: Int,
    val cdsEnd: Int,
    val build: String,
) {
    init {
        require(strand == '+' || strand == '-')
        require(exons.zipWithNext().all { it.first.second < it.second.first }) { "exons must be sorted non-overlapping" }
    }
}

data class PathStep(val step: String, val detail: String)

enum class RegionKind { EXONIC, INTRONIC, UTR5, UTR3, SPLICE_REGION, INTERGENIC, MIXED_EXON_INTRON }

data class ProjectionResult(
    val txId: String,
    val status: ProjectionStatus,
    val failureType: String?,
    val region: RegionKind?,
    val transcriptPos: Int?,
    val refCdna: String?,
    val altCdna: String?,
    val hgvsgLike: String?,
    val path: List<PathStep>,
)

enum class ProjectionStatus { PROJECTED, FAILED }

data class ConsequenceRow(
    val txId: String,
    val txBuild: String,
    val status: String,
    val failureType: String?,
    val region: String?,
    val transcriptPos: Int?,
    val refCdna: String?,
    val altCdna: String?,
    val hgvsgLike: String?,
    val pathJson: String,
)

data class RawVcfRecord(
    val id: Long,
    val batchId: String,
    val source: String,
    val lineNumber: Int,
    val chrom: String,
    val pos: Int,
    val ref: String,
    val altField: String,
    val info: String,
    val rawLine: String,
)

data class BatchReceipt(
    val batchId: String,
    val build: String,
    val receivedAt: String,
    val records: Int,
    val alleles: Int,
    val groups: Int,
    val refMismatches: Int,
    val replay: Boolean,
)

data class Segment(
    val id: Long?,
    val fromBuild: String,
    val toBuild: String,
    val chrom: String,
    val srcStart: Int,
    val srcEnd: Int,
    val strand: Char,
    val dstChrom: String,
    val dstStart: Int,
    val dstEnd: Int,
)

enum class LiftoverStatus { UNIQUE, INTERRUPTED, FLIPPED, MANY_TO_ONE, UNMAPPED }

data class LiftoverResult(
    val status: LiftoverStatus,
    val chrom: String?,
    val pos: Int?,
    val ref: String?,
    val alt: String?,
    val detail: String,
    val candidateCount: Int,
)

data class ImportImpact(
    val fromBuild: String,
    val toBuild: String,
    val segments: Int,
    val newUniqueGroups: Int,
    val flippedGroups: Int,
    val manyToOneGroups: Int,
    val interruptedGroups: Int,
    val unmappedGroups: Int,
    val warnings: List<String>,
)
