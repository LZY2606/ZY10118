package app

/** Failure codes (also documented in README). */
object ProjectionFailures {
    const val UNKNOWN_CHROM = "UNKNOWN_CHROM"
    const val OUT_OF_TRANSCRIPT = "OUT_OF_TRANSCRIPT"
    const val EXON_STRADDLE = "EXON_STRADDLE"
    const val CDS_INCONSISTENT = "CDS_INCONSISTENT"
    const val WRONG_BUILD = "WRONG_BUILD"
}

object TranscriptIo {
    fun parseTsv(content: String): Pair<String, List<Transcript>> {
        var setVersion = "unknown"
        val out = mutableListOf<Transcript>()
        content.lineSequence().forEach { raw ->
            val line = raw.trim()
            if (line.isEmpty()) return@forEach
            if (line.startsWith("#")) {
                Regex("transcriptSetVersion=([^\\s]+)").find(line)?.let { setVersion = it.groupValues[1] }
                return@forEach
            }
            val f = line.split('\t')
            if (f[0] == "txId") return@forEach
            val exons = f[4].split(',').map {
                val (a, b) = it.split('-')
                a.toInt() to b.toInt()
            }
            out += Transcript(f[0], f[1], f[2], f[3].single(), exons, f[5].toInt(), f[6].toInt(), f[7])
        }
        return setVersion to out
    }
}

object TranscriptProjector {

    private val complement = mapOf('A' to 'T', 'T' to 'A', 'C' to 'G', 'G' to 'C')
    fun reverseComplement(s: String): String =
        s.uppercase().reversed().map { complement[it] ?: it }.joinToString("")

    private class TxIndex(tx: Transcript) {
        val segs: List<Triple<Int, Int, Int>>
        val total: Int
        init {
            val list = mutableListOf<Triple<Int, Int, Int>>()
            var cum = 0
            for ((a, b) in tx.exons) {
                list += Triple(a, b, cum)
                cum += b - a + 1
            }
            segs = list
            total = cum
        }
        val txStart get() = segs.first().first
        val txEnd get() = segs.last().second
    }

    /** 1-based cDNA position of an exonic genomic base. */
    private fun g2c(tx: Transcript, idx: TxIndex, gpos: Int): Int {
        val seg = idx.segs.first { gpos in it.first..it.second }
        val forward = gpos - seg.first + 1
        return if (tx.strand == '+') seg.third + forward
        else idx.total - (seg.third + forward) + 1
    }

    private fun clip(tx: Transcript, idx: TxIndex, lo: Int, hi: Int, ea: Int, eb: Int): Pair<Int, Int> {
        val a = maxOf(lo, ea); val b = minOf(hi, eb)
        val ca = g2c(tx, idx, a); val cb = g2c(tx, idx, b)
        return minOf(ca, cb) to maxOf(ca, cb)
    }

    fun project(tx: Transcript, build: String, p: Placement, refGenome: RefGenome): ProjectionResult {
        val path = mutableListOf<PathStep>()
        path += PathStep("input", "placement ${p.chrom}:${p.pos}:${p.ref}>${p.alt} on build $build")
        if (build != tx.build) {
            path += PathStep("fail", "transcript ${tx.txId} is annotated on ${tx.build}, variant is on $build")
            return ProjectionResult(tx.txId, ProjectionStatus.FAILED, ProjectionFailures.WRONG_BUILD,
                null, null, null, null, null, path)
        }
        if (p.chrom != tx.chrom) {
            path += PathStep("fail", "variant contig ${p.chrom} != transcript contig ${tx.chrom}")
            return ProjectionResult(tx.txId, ProjectionStatus.FAILED, ProjectionFailures.UNKNOWN_CHROM,
                null, null, null, null, null, path)
        }
        val idx = TxIndex(tx)
        path += PathStep("exon_index",
            "exons ${tx.exons.joinToString(",") { "${it.first}-${it.second}" }} strand=${tx.strand} cDNA length=${idx.total}")

        val gLo = p.pos
        val gHi = p.pos + p.ref.length - 1
        path += PathStep("affected_interval", "REF spans genomic $gLo..$gHi (1-based closed)")

        if (gHi < idx.txStart || gLo > idx.txEnd) {
            path += PathStep("fail", "interval outside transcript extent ${idx.txStart}..${idx.txEnd}")
            return ProjectionResult(tx.txId, ProjectionStatus.FAILED, ProjectionFailures.OUT_OF_TRANSCRIPT,
                null, null, null, null, null, path)
        }

        data class Hit(val i: Int, val ea: Int, val eb: Int)
        val hits = tx.exons.mapIndexedNotNull { i, (ea, eb) ->
            if (gHi >= ea && gLo <= eb) Hit(i, ea, eb) else null
        }

        if (hits.isEmpty()) {
            var d3 = Int.MAX_VALUE; var d5 = Int.MAX_VALUE
            for ((e1, e2) in tx.exons.zipWithNext()) {
                if (gLo in (e1.second + 1) until e2.first) {
                    d3 = gLo - e1.second
                    d5 = e2.first - gHi
                }
            }
            val near = d3 <= 3 || d5 <= 3
            path += PathStep("classify",
                if (near) "intronic but within 3 bp of a splice junction (upstream exon dist=$d3, downstream dist=$d5)"
                else "no exon overlap; deep intronic")
            path += PathStep("done", "no contiguous cDNA replacement can be assembled for an intronic allele")
            return ProjectionResult(tx.txId, ProjectionStatus.PROJECTED, null,
                if (near) RegionKind.SPLICE_REGION else RegionKind.INTRONIC,
                null, null, null, null, path)
        }

        if (hits.size >= 2) {
            val fe = hits.first(); val se = hits[1]
            path += PathStep("fail",
                "REF $gLo..$gHi spans intron ${fe.eb + 1}..${se.ea - 1} between exons ${fe.i + 1} and ${se.i + 1}")
            return ProjectionResult(tx.txId, ProjectionStatus.FAILED, ProjectionFailures.EXON_STRADDLE,
                null, null, null, null, null, path)
        }

        val h = hits.single()
        // Partial overlap of one exon with an extension into the adjacent
        // intron means the allele removes sequence from both compartments and
        // has no single cDNA representation.
        val crossesIntoIntron = gLo < h.ea || gHi > h.eb
        if (crossesIntoIntron) {
            val where = buildString {
                if (gLo < h.ea) append("${h.ea - gLo} bp upstream intron ")
                if (gHi > h.eb) append("${gHi - h.eb} bp downstream intron")
            }.trim()
            path += PathStep("fail",
                "REF $gLo..$gHi overlaps exon ${h.ea}-${h.eb} but also removes $where; mixed exon/intron allele")
            return ProjectionResult(tx.txId, ProjectionStatus.FAILED, ProjectionFailures.EXON_STRADDLE,
                RegionKind.MIXED_EXON_INTRON, null, null, null, null, path)
        }
        val fullyInside = true
        val (cStart, cEnd) = clip(tx, idx, gLo, gHi, h.ea, h.eb)
        if (!fullyInside) path += PathStep("clip", "partial exon overlap clipped to cDNA $cStart..$cEnd")

        val cdsCdsStart = if (tx.strand == '+') g2c(tx, idx, tx.cdsStart) else g2c(tx, idx, tx.cdsEnd)
        val cdsCdsEnd = if (tx.strand == '+') g2c(tx, idx, tx.cdsEnd) else g2c(tx, idx, tx.cdsStart)
        if (cdsCdsStart > cdsCdsEnd) {
            path += PathStep("fail", "CDS annotation inconsistent with exon structure ($cdsCdsStart > $cdsCdsEnd)")
            return ProjectionResult(tx.txId, ProjectionStatus.FAILED, ProjectionFailures.CDS_INCONSISTENT,
                null, null, null, null, null, path)
        }
        val region = when {
            cEnd < cdsCdsStart || cStart > cdsCdsEnd ->
                if (cStart < cdsCdsStart) RegionKind.UTR5 else RegionKind.UTR3
            else -> RegionKind.EXONIC
        }
        path += PathStep("classify",
            "cDNA $cStart..$cEnd vs CDS cDNA $cdsCdsStart..$cdsCdsEnd -> $region")

        val refWindow = buildString {
            for (g in gLo..gHi) append(refGenome.base(p.chrom, g))
        }
        val refOnStrand = if (tx.strand == '-') reverseComplement(refWindow) else refWindow.uppercase()
        val altOnStrand = if (tx.strand == '-') reverseComplement(p.alt) else p.alt.uppercase()
        path += PathStep("strand",
            if (tx.strand == '-') "negative strand: ${p.ref}>${p.alt} reverse-complemented to $refOnStrand>$altOnStrand"
            else "positive strand: allele used as written ($refOnStrand>$altOnStrand)")

        val hgvsg = "${tx.txId}:c.${describeHgvs(cStart, cEnd, refOnStrand, altOnStrand)}"
        path += PathStep("hgvsg", "assembled transcript-relative notation $hgvsg")
        path += PathStep("done", "projection complete with rules $ANNOTATION_RULES_VERSION / transcript set $TRANSCRIPT_SET_VERSION")
        return ProjectionResult(tx.txId, ProjectionStatus.PROJECTED, null, region,
            cStart, refOnStrand, altOnStrand, hgvsg, path)
    }

    private fun describeHgvs(cStart: Int, cEnd: Int, refS: String, altS: String): String = when {
        refS.length == 1 && altS.length == 1 -> "$cStart$refS>$altS"
        altS.length > refS.length -> "${cStart}_${cEnd}ins${altS.drop(1)}"
        refS.length > altS.length -> "${cStart}_${cEnd}del${refS.drop(1)}"
        else -> "${cStart}_${cEnd}del${refS}ins$altS"
    }
}
