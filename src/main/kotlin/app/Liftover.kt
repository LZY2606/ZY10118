package app

/**
 * Minimal chain fixture for coordinate mapping between reference builds.
 * TSV columns (all 1-based closed):
 *   fromBuild toBuild chrom srcStart srcEnd strand dstChrom dstStart dstEnd
 */
object ChainIo {
    fun parse(content: String): Pair<Pair<String, String>, List<Segment>> {
        var from = ""; var to = ""
        val segs = mutableListOf<Segment>()
        content.lineSequence().forEach { raw ->
            val line = raw.trim()
            if (line.isEmpty()) return@forEach
            if (line.startsWith("#")) {
                Regex("fromBuild=(\\S+)").find(line)?.let { from = it.groupValues[1] }
                Regex("toBuild=(\\S+)").find(line)?.let { to = it.groupValues[1] }
                return@forEach
            }
            val f = line.split('\t')
            require(f.size == 9) { "chain row must have 9 columns, got ${f.size}: $line" }
            segs += Segment(
                id = null, fromBuild = f[0], toBuild = f[1], chrom = f[2],
                srcStart = f[3].toInt(), srcEnd = f[4].toInt(), strand = f[5].single(),
                dstChrom = f[6], dstStart = f[7].toInt(), dstEnd = f[8].toInt(),
            )
        }
        require(segs.isNotEmpty()) { "chain file contains no segments" }
        require(segs.all { it.srcEnd - it.srcStart == it.dstEnd - it.dstStart }) {
            "chain segment source/destination lengths must agree"
        }
        return (from to to) to segs
    }
}

object LiftoverEngine {

    fun map(p: Placement, segs: List<Segment>): LiftoverResult {
        val spanEnd = p.pos + p.ref.length - 1
        val covering = segs.filter {
            it.chrom == p.chrom && p.pos >= it.srcStart && spanEnd <= it.srcEnd
        }
        if (covering.isEmpty()) {
            // Distinguish "no segment anywhere near" (UNMAPPED) from "segments
            // touch the interval but the chain is broken across its span"
            // (INTERRUPTED).
            val touching = segs.filter {
                it.chrom == p.chrom && it.srcEnd >= p.pos && it.srcStart <= spanEnd
            }
            if (touching.isNotEmpty()) {
                val names = touching.joinToString(",") { "${it.srcStart}-${it.srcEnd}" }
                return LiftoverResult(LiftoverStatus.INTERRUPTED, null, null, null, null,
                    "interval ${p.pos}..$spanEnd is only partially covered by chain segment(s) $names; chain breaks inside the allele span",
                    touching.size)
            }
            val near = segs.filter { it.chrom == p.chrom }.minByOrNull {
                minOf(kotlin.math.abs(p.pos - it.srcStart), kotlin.math.abs(p.pos - it.srcEnd))
            }
            return LiftoverResult(
                LiftoverStatus.UNMAPPED, null, null, null, null,
                if (near != null) "no chain segment covers ${p.chrom}:${p.pos}..$spanEnd; nearest source interval ${near.srcStart}-${near.srcEnd}"
                else "no chain segments for contig ${p.chrom}", 0,
            )
        }
        if (covering.size > 1) {
            // The interval is claimed by more than one source segment. If the
            // segments map the anchor to different destination coordinates (or
            // disagree in direction), the source position has no single image:
            // MANY_TO_ONE. Identical images are treated as redundant coverage.
            data class Image(val strand: Char, val dstChrom: String, val dstPos: Int)
            fun imageOf(o: Segment): Image {
                val off = p.pos - o.srcStart
                val dp = if (o.strand == '+') o.dstStart + off else o.dstEnd - off
                return Image(o.strand, o.dstChrom, dp)
            }
            val images = covering.map(::imageOf).distinct()
            val flippedConflict = images.map { it.strand }.distinct().size > 1
            return if (images.size == 1 && !flippedConflict) {
                LiftoverResult(LiftoverStatus.INTERRUPTED, null, null, null, null,
                    "${covering.size} segments claim ${p.pos}..$spanEnd but agree on the image; inspect chain",
                    covering.size)
            } else {
                val detail = images.joinToString(",") { "${it.strand} ${it.dstChrom}:${it.dstPos}" }
                LiftoverResult(LiftoverStatus.MANY_TO_ONE, null, null, null, null,
                    "${covering.size} chain segments map ${p.chrom}:${p.pos} to distinct destinations ($detail)",
                    covering.size)
            }
        }
        val seg = covering.single()
        val flipped = seg.strand == '-'
        val offset = p.pos - seg.srcStart
        val dstPos = if (!flipped) seg.dstStart + offset else seg.dstEnd - offset
        val refT = if (flipped) TranscriptProjector.reverseComplement(p.ref) else p.ref
        val altT = if (flipped) TranscriptProjector.reverseComplement(p.alt) else p.alt

        if (spanEnd == seg.srcEnd && p.ref.length > 1) {
            return LiftoverResult(
                LiftoverStatus.INTERRUPTED, seg.dstChrom,
                if (!flipped) seg.dstEnd else seg.dstStart, refT, altT,
                "allele REF ends exactly at chain segment boundary ${seg.srcEnd}; no continuation segment, mapping interrupted", 1,
            )
        }

        // MANY_TO_ONE: destination intervals of two distinct source segments
        // overlap while neither interval strictly contains the other. (A small
        // interval tiled inside a larger one does not by itself make coordinates
        // ambiguous for variants inside the large interval.)
        fun overlaps(a: Segment, b: Segment) =
            a.dstChrom == b.dstChrom && a.dstStart <= b.dstEnd && b.dstStart <= a.dstEnd
        fun contains(a: Segment, b: Segment) =
            a.dstChrom == b.dstChrom && a.dstStart <= b.dstStart && b.dstEnd <= a.dstEnd
        val collisions = segs.filter { o ->
            o !== seg && overlaps(seg, o) && !contains(seg, o) && !contains(o, seg)
        }
        if (collisions.isNotEmpty()) {
            val names = collisions.joinToString(",") { "${it.chrom}:${it.srcStart}-${it.srcEnd}" }
            return LiftoverResult(LiftoverStatus.MANY_TO_ONE, seg.dstChrom, dstPos, refT, altT,
                "destination interval ${seg.dstStart}-${seg.dstEnd} collides with destination(s) from $names; coordinate not unique",
                collisions.size + 1)
        }
        if (flipped) {
            return LiftoverResult(LiftoverStatus.FLIPPED, seg.dstChrom, dstPos, refT, altT,
                "source ${seg.srcStart}-${seg.srcEnd}(+) maps reversed to ${seg.dstChrom}:${seg.dstStart}-${seg.dstEnd}(-)", 1)
        }
        return LiftoverResult(LiftoverStatus.UNIQUE, seg.dstChrom, dstPos, refT, altT,
            "unique forward-strand mapping: ${seg.srcStart}-${seg.srcEnd} -> ${seg.dstStart}-${seg.dstEnd}", 1)
    }
}
