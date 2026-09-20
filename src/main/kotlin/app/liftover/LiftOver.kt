package app.liftover

import app.model.LiftResult
import app.model.MapSegment
import app.model.Normalized
import java.io.File

/** Lift status values (documented in README):
 *  MAPPED          unique same-strand placement
 *  STRAND_FLIPPED  unique placement but destination chain is on the opposite strand
 *  GAP             the allele span falls in an unmapped assembly gap
 *  COLLISION       more than one chain claims the source region (many-to-one)
 *  OUT_OF_MAP      no chain covers the source region at all
 */
class LiftOver(
    val sourceBuild: String,
    val destBuild: String,
    val segments: List<MapSegment>,
) {
    /** Segments covering [start, end] (1-based inclusive) on the source contig. */
    private fun covering(contig: String, start: Int, end: Int): List<Pair<Int, MapSegment>> =
        segments.withIndex()
            .filter { (_, s) -> s.srcContig == contig && s.srcStart <= end && start <= s.srcEnd }
            .map { it.index + 1 to it.value }

    fun lift(contig: String, allele: Normalized, destReference: Map<String, String>): LiftResult {
        val start = allele.pos
        val end = allele.pos + allele.ref.length - 1
        val hits = covering(contig, start, end)
        if (hits.isEmpty()) {
            return LiftResult("OUT_OF_MAP", null, null, null, null, emptyList(),
                "no chain covers $contig:$start-$end")
        }
        // A true many-to-one collision requires two segments whose SOURCE ranges
        // overlap. Adjacent (non-overlapping) segments merely share a boundary;
        // an allele spanning that break crosses an assembly gap.
        val overlap = hits.any { (i1, s1) ->
            hits.any { (i2, s2) -> i1 != i2 && s1.srcStart <= s2.srcEnd && s2.srcStart <= s1.srcEnd &&
                // overlapping beyond a single touching base means duplicated source
                !(s1.srcEnd + 1 == s2.srcStart || s2.srcEnd + 1 == s1.srcStart) }
        }
        if (overlap) {
            return LiftResult(
                "COLLISION", null, null, null, null,
                hits.map { it.first },
                "${hits.size} chains claim overlapping source at $contig:$start (many-to-one); lift not unique"
            )
        }
        if (hits.size > 1) {
            return LiftResult("GAP", null, null, null, null, hits.map { it.first },
                "allele spans the assembly gap between segments " +
                    hits.joinToString("#") { it.first.toString() } + " at $contig:$start-$end")
        }
        val (segIndex, seg) = hits.single()
        // Span must lie fully inside one segment, otherwise it crosses a gap.
        if (start < seg.srcStart || end > seg.srcEnd) {
            return LiftResult("GAP", null, null, null, null, listOf(segIndex),
                "allele spans an assembly gap at the edge of segment #$segIndex")
        }

        val offset = if (seg.dstStrand == "+") {
            start - seg.srcStart + seg.dstStart
        } else {
            seg.dstEnd - (start - seg.srcStart)
        }
        val destContig = seg.dstContig
        val destSeq = destReference[destContig]
            ?: return LiftResult("OUT_OF_MAP", null, null, null, null, listOf(segIndex),
                "destination contig $destContig is absent from build $destBuild")

        val destRef: String
        val destAlt: String
        if (seg.dstStrand == "+") {
            destRef = allele.ref
            destAlt = allele.alt
        } else {
            destRef = reverseComplement(allele.ref)
            destAlt = reverseComplement(allele.alt)
        }
        if (offset - 1 + destRef.length > destSeq.length) {
            return LiftResult("GAP", null, null, null, null, listOf(segIndex),
                "lifted allele runs past the destination contig end")
        }
        val observed = destSeq.substring(offset - 1, offset - 1 + destRef.length)
        if (!observed.equals(destRef, ignoreCase = true)) {
            return LiftResult("GAP", null, null, null, null, listOf(segIndex),
                "destination bases $observed do not match lifted REF $destRef")
        }
        val status = if (seg.dstStrand == "-") "STRAND_FLIPPED" else "MAPPED"
        val note = if (status == "STRAND_FLIPPED")
            "segment #$segIndex maps to the opposite strand; alleles reverse-complemented"
        else "mapped via segment #$segIndex with offset ${seg.dstStart - seg.srcStart}"
        return LiftResult(status, destContig, offset, destRef, destAlt, listOf(segIndex), note)
    }

    companion object {
        fun reverseComplement(s: String): String {
            val comp = mapOf('A' to 'T', 'T' to 'A', 'C' to 'G', 'G' to 'C')
            return s.reversed().map { comp[it] ?: it }.joinToString("")
        }

        fun parse(file: File, sourceBuild: String, destBuild: String): LiftOver {
            val rows = file.readLines().filter { it.isNotBlank() }
            require(rows.size >= 2) { "chain TSV ${file.name} has no segments" }
            val segments = rows.drop(1).map { line ->
                val c = line.split("\t")
                require(c.size == 7) { "chain row needs 7 columns: $line" }
                MapSegment(
                    srcContig = c[0],
                    srcStart = c[1].toInt(), srcEnd = c[2].toInt(),
                    dstContig = c[3],
                    dstStart = c[4].toInt(), dstEnd = c[5].toInt(),
                    dstStrand = c[6],
                )
            }
            return LiftOver(sourceBuild, destBuild, segments)
        }
    }
}
