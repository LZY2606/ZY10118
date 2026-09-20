package app.transcript

import app.model.ANNOTATION_RULES_VERSION
import app.model.Normalized
import app.model.Projection
import app.model.Transcript
import kotlin.math.abs

/** Failure codes used when a consequence cannot be uniquely assigned. They are
 *  part of the public contract documented in README.md. */
object FailureType {
    const val CONTIG_NOT_FOUND = "CONTIG_NOT_FOUND"
    const val OUTSIDE_TRANSCRIPT = "OUTSIDE_TRANSCRIPT"
    const val SPANS_INTRON = "SPANS_INTRON"
    const val REFERENCE_MISMATCH = "REFERENCE_MISMATCH"
    const val NON_UNIQUE_LEFT_ALIGN = "NON_UNIQUE_LEFT_ALIGN"
    const val UNKNOWN_ALLELE = "UNKNOWN_ALLELE"
}

/**
 * Projects a canonical allele onto one transcript.
 *
 * Allele span model after normalization (1-based):
 *   - the anchor base lives at pos and is identical in REF/ALT;
 *   - removed reference bases occupy pos+1 .. pos+len(ref)-1;
 *   - inserted bases (len(alt)-1 of them) are attached after the anchor.
 *
 * A substitution (len 1/1) changes the base exactly at pos.
 */
class Projector(
    private val buildName: String,
    private val reference: Map<String, String>,
) {
    private val codonTable = mapOf(
        "TTT" to "F", "TTC" to "F", "TTA" to "L", "TTG" to "L",
        "CTT" to "L", "CTC" to "L", "CTA" to "L", "CTG" to "L",
        "ATT" to "I", "ATC" to "I", "ATA" to "I", "ATG" to "M",
        "GTT" to "V", "GTC" to "V", "GTA" to "V", "GTG" to "V",
        "TCT" to "S", "TCC" to "S", "TCA" to "S", "TCG" to "S",
        "CCT" to "P", "CCC" to "P", "CCA" to "P", "CCG" to "P",
        "ACT" to "T", "ACC" to "T", "ACA" to "T", "ACG" to "T",
        "GCT" to "A", "GCC" to "A", "GCA" to "A", "GCG" to "A",
        "TAT" to "Y", "TAC" to "Y", "TAA" to "*", "TAG" to "*",
        "CAT" to "H", "CAC" to "H", "CAA" to "Q", "CAG" to "Q",
        "AAT" to "N", "AAC" to "N", "AAA" to "K", "AAG" to "K",
        "GAT" to "D", "GAC" to "D", "GAA" to "E", "GAG" to "E",
        "TGT" to "C", "TGC" to "C", "TGA" to "*", "TGG" to "W",
        "CGT" to "R", "CGC" to "R", "CGA" to "R", "CGG" to "R",
        "AGT" to "S", "AGC" to "S", "AGA" to "R", "AGG" to "R",
        "GGT" to "G", "GGC" to "G", "GGA" to "G", "GGG" to "G",
    )

    fun project(tx: Transcript, allele: Normalized, refMatches: Boolean): Projection {
        val path = ArrayList<String>()
        val seq = reference[tx.contig]
        if (seq == null) return fail(tx, path, FailureType.CONTIG_NOT_FOUND,
            "contig ${tx.contig} is absent from build $buildName")
        if (!refMatches) return fail(tx, path, FailureType.REFERENCE_MISMATCH,
            "REF did not match build $buildName; projection withheld")
        if (!allele.uniquePlacement) return fail(tx, path,
            FailureType.NON_UNIQUE_LEFT_ALIGN,
            "allele has ${allele.equivalentAnchors.size} equivalent physical anchors")
        if ((allele.ref + allele.alt).any { it !in "ACGT" })
            return fail(tx, path, FailureType.UNKNOWN_ALLELE, "allele contains ambiguity codes")

        path.add("allele canonical at ${tx.contig}:${allele.pos} ${allele.ref}>${allele.alt}")

        val isSub = allele.ref.length == 1 && allele.alt.length == 1
        // Positions whose reference bases are removed (empty for a substitution).
        val removedPositions = (allele.pos + 1 until allele.pos + allele.ref.length).toList()
        // For substitutions the changed base is the anchor itself.
        val touchedPositions: List<Int> =
            if (isSub) listOf(allele.pos) else removedPositions + listOf(allele.pos)

        val spanLo = allele.pos
        val spanHi = allele.pos + allele.ref.length - 1
        // The transcript extent runs from the first exon start to the last
        // exon end (introns included); only points beyond that are intergenic.
        val txStart = tx.exons.minOf { it.first }
        val txEnd = tx.exons.maxOf { it.second }
        if (spanHi < txStart || spanLo > txEnd) {
            return fail(tx, path, FailureType.OUTSIDE_TRANSCRIPT,
                "allele span $spanLo-$spanHi is outside transcript ${tx.id} ($txStart-$txEnd)")
        }

        val changedPositions = if (isSub) listOf(allele.pos) else removedPositions
        val exonsHit = changedPositions.filter { tx.inExon(it) }
            .mapNotNull { exonIndex(tx, it) }.distinct().sorted()

        // Pure intronic: the anchor sits inside an intron and nothing exonic changed.
        val anchorExonic = tx.inExon(allele.pos)
        if (!isSub && changedPositions.none { tx.inExon(it) } && !anchorExonic) {
            path.add("anchor and all ${removedPositions.size} removed base(s) are intronic")
            return consequence(tx, path, "INTRON_VARIANT", exonsHit)
        }
        if (isSub && !tx.inExon(allele.pos)) {
            path.add("substitution base ${allele.pos} is intronic")
            return consequence(tx, path, "INTRON_VARIANT", emptyList())
        }
        // An indel whose removed bases straddle an exon edge cannot be assigned
        // a single coding outcome.
        if (!isSub && removedPositions.isNotEmpty() &&
            removedPositions.any { tx.inExon(it) } &&
            removedPositions.any { !tx.inExon(it) }
        ) {
            return fail(tx, path, FailureType.SPANS_INTRON,
                "indel removed bases cross the exon/intron boundary at ${tx.id}")
        }
        // Insertion anchored exactly on an exon/intron boundary.
        if (!isSub && removedPositions.isEmpty() &&
            tx.exons.any { allele.pos == it.first || allele.pos == it.second }
        ) {
            return fail(tx, path, FailureType.SPANS_INTRON,
                "insertion anchor sits on an exon boundary; splicing impact ambiguous")
        }

        val codingChanged = changedPositions.filter { tx.inCoding(it) }
        if (codingChanged.isEmpty()) {
            val label = utrLabel(tx, touchedPositions)
            path.add("no coding base changed; classified as $label")
            return consequence(tx, path, label, exonsHit)
        }

        val cdsCoordinates = codingChanged.mapNotNull { cdsCoordinate(tx, it) }.sorted()
        val cdsFirst = cdsCoordinates.first()
        val cdsLast = cdsCoordinates.last()
        path.add("coding change covers CDS $cdsFirst-$cdsLast (strand ${tx.strand})")

        if (isSub) {
            val g = allele.pos
            val cds = cdsCoordinate(tx, g)!!
            val codonIndex0 = (cds - 1) / 3
            val refCodon = readCodon(tx, seq, codonIndex0)
            val altCodonChars = refCodon.toCharArray()
            altCodonChars[(cds - 1) % 3] = allele.alt.single()
            val altCodon = String(altCodonChars)
            val refAa = codonTable[refCodon] ?: "?"
            val altAa = codonTable[altCodon] ?: "?"
            path.add("CDS $cds belongs to codon ${codonIndex0 + 1} offset ${(cds - 1) % 3 + 1}")
            path.add("codon $refCodon($refAa) -> $altCodon($altAa)")
            val name = when {
                codonIndex0 == 0 && (cds in 1..3) && refAa == "M" && altAa != "M" -> "START_LOST"
                refAa == "*" && altAa != "*" -> "STOP_LOST"
                refAa != "*" && altAa == "*" -> "STOP_GAINED"
                refAa == altAa -> "SYNONYMOUS_VARIANT"
                else -> "MISSENSE_VARIANT"
            }
            return Projection(
                transcriptId = tx.id, build = buildName,
                rulesVersion = ANNOTATION_RULES_VERSION, consequence = name, failure = null,
                cdsStart = cds, cdsEnd = cds,
                proteinStart = codonIndex0 + 1, proteinEnd = codonIndex0 + 1,
                refCodon = refCodon, altCodon = altCodon, refAa = refAa, altAa = altAa,
                path = path, affectedExons = exonsHit,
            )
        }

        val removed = allele.ref.length - 1
        val inserted = allele.alt.length - 1
        val net = inserted - removed
        path.add("removed $removed reference base(s), inserted $inserted (net $net)")
        val proteinStart = (cdsFirst - 1) / 3 + 1
        val name = if (abs(net) % 3 != 0) "FRAMESHIFT_VARIANT" else "INFRAME_INDEL"
        return Projection(
            transcriptId = tx.id, build = buildName,
            rulesVersion = ANNOTATION_RULES_VERSION, consequence = name, failure = null,
            cdsStart = cdsFirst, cdsEnd = cdsLast,
            proteinStart = proteinStart, proteinEnd = proteinStart,
            refCodon = null, altCodon = null, refAa = null, altAa = null,
            path = path, affectedExons = exonsHit,
        )
    }

    private fun exonIndex(tx: Transcript, g: Int): Int? =
        tx.exons.indexOfFirst { it.first <= g && g <= it.second }
            .let { if (it < 0) null else it + 1 }

    private fun utrLabel(tx: Transcript, positions: List<Int>): String {
        val fivePrime = if (tx.strand == "+")
            positions.any { it < tx.cdsStart } else positions.any { it > tx.cdsEnd }
        return if (fivePrime) "FIVE_PRIME_UTR_VARIANT" else "THREE_PRIME_UTR_VARIANT"
    }

    fun cdsCoordinate(tx: Transcript, g: Int): Int? {
        if (!tx.inCoding(g)) return null
        return if (tx.strand == "+") {
            var count = 0
            for ((s, e) in tx.codingPieces) if (g >= s) count += minOf(g, e) - s + 1
            count
        } else {
            var count = 0
            for ((s, e) in tx.codingPieces.asReversed()) if (e >= g) count += e - maxOf(g, s) + 1
            count
        }
    }

    private fun genomicForCds(tx: Transcript, cds1: Int): Int {
        var remaining = cds1
        val pieces = if (tx.strand == "+") tx.codingPieces else tx.codingPieces.asReversed()
        for ((s, e) in pieces) {
            val len = e - s + 1
            if (remaining <= len) {
                return if (tx.strand == "+") s + remaining - 1 else e - remaining + 1
            }
            remaining -= len
        }
        throw IllegalStateException("CDS coordinate $cds1 outside transcript ${tx.id}")
    }

    private fun readCodon(tx: Transcript, seq: String, codonIndex0: Int): String {
        val sb = StringBuilder()
        for (k in 0..2) {
            val g = genomicForCds(tx, codonIndex0 * 3 + k + 1)
            val b = seq[g - 1].uppercaseChar()
            sb.append(if (tx.strand == "-") complement(b) else b)
        }
        return sb.toString()
    }

    private fun complement(b: Char): Char = when (b) {
        'A' -> 'T'; 'T' -> 'A'; 'C' -> 'G'; 'G' -> 'C'; else -> b
    }

    private fun fail(tx: Transcript, path: MutableList<String>, code: String, note: String) =
        Projection(
            transcriptId = tx.id, build = buildName,
            rulesVersion = ANNOTATION_RULES_VERSION,
            consequence = "PROJECTION_FAILED", failure = code,
            cdsStart = null, cdsEnd = null, proteinStart = null, proteinEnd = null,
            refCodon = null, altCodon = null, refAa = null, altAa = null,
            path = path + "failure: $note", affectedExons = emptyList(),
        )

    private fun consequence(tx: Transcript, path: MutableList<String>, name: String,
                            exons: List<Int>) = Projection(
        transcriptId = tx.id, build = buildName,
        rulesVersion = ANNOTATION_RULES_VERSION, consequence = name, failure = null,
        cdsStart = null, cdsEnd = null, proteinStart = null, proteinEnd = null,
        refCodon = null, altCodon = null, refAa = null, altAa = null,
        path = path, affectedExons = exons,
    )
}
