package app

import org.junit.jupiter.api.Assertions.*
import org.junit.jupiter.api.Test

class NormalizerTest {
    /**
     * Reference layout (1-based):
     *  699-708: C A A A A A A A A C  (8-A run, C anchors)
     *  799-808: C G T G T G T G T C  ((GT)4 with C anchors)
     *  899-909: C C A T C A T C A T  ((CAT)3 with C anchors, pos 909 T)
     *  100.. : ACGTACGT
     */
    private val seq = run {
        val sb = StringBuilder("N".repeat(2000))
        fun put(pos1: Int, s: String) {
            for ((i, c) in s.withIndex()) sb.setCharAt(pos1 - 1 + i, c)
        }
        put(699, "CAAAAAAAAC")
        put(799, "CGTGTGTGTC")
        put(899, "CCATCATCAT")
        put(100, "ACGTACGT")
        sb.toString()
    }
    private val genome = InMemoryRefGenome("hg19", seq)

    @Test
    fun `left-aligns insertion in homopolymer to earliest anchor position`() {
        // AA -> AAA entered at 700 (inside the A run) slides left until anchor C@699
        val n = Normalizer.normalize(Variant("chr1", 700, "AA", "AAA"), genome)
        assertTrue(n.refMatches)
        val c = n.canonical!!
        assertEquals(699, c.pos)
        assertEquals("C", c.ref)
        assertEquals("CA", c.alt)
        // A^8 between two C anchors: 8 internal gaps plus both C anchors => 10 minimal placements
        assertEquals(10, n.placements.size)
        val positions = n.placements.map { it.pos }
        assertEquals(positions.sorted(), positions)
    }

    @Test
    fun `same change entered at another run slot gives identical canonical placement`() {
        // At 702 ref bases are A A: AA -> AAA is the same +A event
        val n = Normalizer.normalize(Variant("chr1", 702, "AA", "AAA"), genome)
        val c = n.canonical!!
        assertEquals(699, c.pos)
        assertEquals("C", c.ref)
        assertEquals("CA", c.alt)
    }

    @Test
    fun `tandem repeat GT keeps all equivalent placements`() {
        // +GT entered at 800 with GTG -> GTGTG
        val n = Normalizer.normalize(Variant("chr1", 800, "GTG", "GTGTG"), genome)
        assertTrue(n.placements.size >= 2)
        val c = n.canonical!!
        assertEquals(799, c.pos)
        assertEquals("C", c.ref)
        assertEquals("CGT", c.alt)
        val keys = n.placements.map { "${it.pos}:${it.ref}>${it.alt}" }.toSet()
        assertEquals(n.placements.size, keys.size)
        assertTrue(n.steps.any { it.kind == "repeat_ambiguity" })
    }

    @Test
    fun `tandem deletion left aligns and enumerates slots`() {
        // CATC -> C at 900 (delete CAT inside the (CAT)3 run)
        val n = Normalizer.normalize(Variant("chr1", 900, "CATC", "C"), genome)
        val c = n.canonical!!
        assertEquals(899, c.pos)
        assertEquals("CCAT", c.ref)
        assertEquals("C", c.alt)
        assertTrue(n.placements.size >= 2)
    }

    @Test
    fun `snv is a single placement`() {
        val n = Normalizer.normalize(Variant("chr1", 100, "A", "T"), genome)
        assertEquals(1, n.placements.size)
        assertEquals(100, n.canonical!!.pos)
    }

    @Test
    fun `ref mismatch is reported with expected reference`() {
        val n = Normalizer.normalize(Variant("chr1", 100, "T", "A"), genome)
        assertFalse(n.refMatches)
        assertEquals("A", n.expectedRef)
        assertNull(n.canonical)
    }

    @Test
    fun `out of bounds reference is mismatch not crash`() {
        val n = Normalizer.normalize(Variant("chr1", 100000, "A", "T"), genome)
        assertFalse(n.refMatches)
        assertNull(n.expectedRef)
    }

    @Test
    fun `minimal representation strips shared anchors`() {
        // 100 ACGTACGT: ACGTA -> ACGTT shares prefix ACGT
        val n = Normalizer.normalize(Variant("chr1", 100, "ACGTA", "ACGTT"), genome)
        val c = n.canonical!!
        assertEquals(104, c.pos)
        assertEquals("A", c.ref)
        assertEquals("T", c.alt)
    }
}
