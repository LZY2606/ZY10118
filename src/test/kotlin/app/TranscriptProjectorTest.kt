package app

import org.junit.jupiter.api.Assertions.*
import org.junit.jupiter.api.Test

class TranscriptProjectorTest {
    private val plusTx = Transcript("TX_P", "G", "chr1", '+',
        listOf(201 to 289, 370 to 459), 209, 456, "hg19")
    private val minusTx = Transcript("TX_N", "G2", "chr1", '-',
        listOf(1501 to 1580, 1661 to 1740), 1512, 1732, "hg19")

    // reference covering the transcript regions deterministically
    private fun genomeFor(): InMemoryRefGenome {
        val sb = StringBuilder("N".repeat(2000))
        fun put(pos1: Int, s: String) { for ((i, c) in s.withIndex()) sb.setCharAt(pos1 - 1 + i, c) }
        // exon1 bases are A-except marker, exon2 bases T-ish; intron 290..369 G
        for (p in 201..289) sb.setCharAt(p - 1, 'A')
        for (p in 370..459) sb.setCharAt(p - 1, 'C')
        for (p in 290..369) sb.setCharAt(p - 1, 'G')
        for (p in 1501..1740) sb.setCharAt(p - 1, 'T')
        put(220, "A")  // exon1 SNV anchor
        return InMemoryRefGenome("hg19", sb.toString())
    }
    private val genome = genomeFor()

    @Test
    fun chromosomeStart_snvOnPlusStrand() {
        val r = TranscriptProjector.project(plusTx, "hg19", Placement("chr1", 220, "A", "T", 0), genome)
        assertEquals(ProjectionStatus.PROJECTED, r.status)
        assertEquals(RegionKind.EXONIC, r.region)
        assertEquals(20, r.transcriptPos)
        assertEquals("A", r.refCdna); assertEquals("T", r.altCdna)
        assertTrue(r.hgvsgLike!!.contains(":c.20A>T"))
        assertTrue(r.path.any { it.detail.contains("positive strand") })
    }

    @Test
    fun insertionAndDeleteInsideExonProjectToCdna() {
        // insertion at exon1 base 220: A -> AA
        val ins = TranscriptProjector.project(plusTx, "hg19", Placement("chr1", 220, "A", "AA", 0), genome)
        assertEquals(ProjectionStatus.PROJECTED, ins.status)
        assertEquals("A", ins.refCdna); assertEquals("AA", ins.altCdna)
        // deletion 220..222 AAA -> A
        val del = TranscriptProjector.project(plusTx, "hg19", Placement("chr1", 220, "AAA", "A", 0), genome)
        assertEquals(ProjectionStatus.PROJECTED, del.status)
        assertTrue(del.hgvsgLike!!.contains("del"))
    }

    @Test
    fun deletionAcrossExonBoundaryIsExonStraddle() {
        // 285..294 covers last 5 exon1 bases plus 5 intron bases
        val r = TranscriptProjector.project(plusTx, "hg19", Placement("chr1", 285, "A".repeat(5) + "G".repeat(5), "AAAAA", 0), genome)
        assertEquals(ProjectionStatus.FAILED, r.status)
        assertEquals(ProjectionFailures.EXON_STRADDLE, r.failureType)
    }

    @Test
    fun deepIntronicVsSpliceRegion() {
        val deep = TranscriptProjector.project(plusTx, "hg19", Placement("chr1", 330, "G", "T", 0), genome)
        assertEquals(RegionKind.INTRONIC, deep.region)
        val near = TranscriptProjector.project(plusTx, "hg19", Placement("chr1", 291, "G", "T", 0), genome)
        assertEquals(RegionKind.SPLICE_REGION, near.region)
    }

    @Test
    fun negativeStrandReverseComplementsAllele() {
        // genomic T->A is reverse-complemented to cDNA A->T on the minus strand
        val r = TranscriptProjector.project(minusTx, "hg19", Placement("chr1", 1550, "T", "A", 0), genome)
        assertEquals(ProjectionStatus.PROJECTED, r.status)
        assertEquals(RegionKind.EXONIC, r.region)
        assertEquals("A", r.refCdna)
        assertEquals("T", r.altCdna)
        assertTrue(r.path.any { it.detail.contains("reverse-complemented") })
    }

    @Test
    fun outOfTranscriptAndWrongBuildFailures() {
        val out = TranscriptProjector.project(plusTx, "hg19", Placement("chr1", 1000, "G", "T", 0), genome)
        assertEquals(ProjectionFailures.OUT_OF_TRANSCRIPT, out.failureType)
        val wb = TranscriptProjector.project(plusTx, "hg38", Placement("chr1", 220, "A", "T", 0), genome)
        assertEquals(ProjectionFailures.WRONG_BUILD, wb.failureType)
    }
}
