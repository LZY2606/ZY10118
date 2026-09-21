package app

import org.junit.jupiter.api.Assertions.*
import org.junit.jupiter.api.Test

class LiftoverTest {
    private val segs = listOf(
        Segment(1, "hg19", "hg38", "chr1", 100, 1999, '+', "chr1", 1000, 2899),
        Segment(2, "hg19", "hg38", "chr1", 2000, 2099, '-', "chr1", 4000, 4099),
        Segment(3, "hg19", "hg38", "chr1", 2010, 2099, '-', "chr1", 5050, 5139),
        Segment(4, "hg19", "hg38", "chr1", 2200, 2240, '+', "chr1", 6000, 6040),
    )
    private fun place(pos: Int, ref: String = "A", alt: String = "T") = Placement("chr1", pos, ref, alt, 0)

    @Test
    fun uniqueForwardMapping() {
        val r = LiftoverEngine.map(place(100), segs)
        assertEquals(LiftoverStatus.UNIQUE, r.status)
        assertEquals(1000, r.pos)
    }

    @Test
    fun flippedSegmentReversesCoordinate() {
        val r = LiftoverEngine.map(place(2000), segs.filter { it.id == 2L })
        assertEquals(LiftoverStatus.FLIPPED, r.status)
        assertEquals(4099, r.pos)
    }

    @Test
    fun overlappingSegmentsWithDifferentDestinationsAreManyToOne() {
        val r = LiftoverEngine.map(place(2050), segs)
        assertEquals(LiftoverStatus.MANY_TO_ONE, r.status)
        assertTrue(r.detail.contains("distinct destinations"))
    }

    @Test
    fun gapInsideAlleleSpanIsInterrupted() {
        // 2238..2241 deletion; seg4 ends at 2240, no continuation until 2250
        val r = LiftoverEngine.map(place(2238, "GAGT", "G"), segs)
        assertEquals(LiftoverStatus.INTERRUPTED, r.status)
    }

    @Test
    fun uncoveredPositionIsUnmapped() {
        val r = LiftoverEngine.map(place(2500), segs)
        assertEquals(LiftoverStatus.UNMAPPED, r.status)
        assertEquals(0, r.candidateCount)
    }

    @Test
    fun parserRejectsLengthMismatch() {
        val bad = "# fromBuild=hg19 toBuild=hg38 segments=1\nhg19\thg38\tchr1\t1\t10\t+\tchr1\t100\t110\n"
        assertThrows(IllegalArgumentException::class.java) { ChainIo.parse(bad) }
    }
}
