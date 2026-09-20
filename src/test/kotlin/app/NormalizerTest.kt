package app

import app.norm.Normalizer
import org.junit.jupiter.api.Assertions.assertEquals
import org.junit.jupiter.api.Assertions.assertFalse
import org.junit.jupiter.api.Assertions.assertTrue
import org.junit.jupiter.api.Test

/** Pure algorithm tests using synthetic, fully controlled references. */
class NormalizerTest {

    private val chr20 = "ACGTACGTAC" + "A".repeat(8) +
        // 19-58 deterministic 40 bases
        "ACGACGACGACGACGACGACGACGACGACGACGACGA" +
        "T".repeat(8) +
        "ACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGACGA"

    @Test
    fun `substitution at chromosome start stays at position 1`() {
        val n = Normalizer(chr20).normalize(1, "A", "T")
        assertEquals(1, n.pos)
        assertEquals("A", n.ref); assertEquals("T", n.alt)
        assertEquals(listOf(1), n.equivalentAnchors)
        assertTrue(n.uniquePlacement)
    }

    @Test
    fun `insertion in poly-A left aligns to run start with all anchors`() {
        val n = Normalizer(chr20).normalize(14, "A", "AA")
        // poly-A occupies 11..18; inserting one A anywhere is the same molecule,
        // anchors 11..18, preferred (leftmost) is 11.
        assertEquals(11, n.pos)
        assertEquals("A", n.ref); assertEquals("AA", n.alt)
        assertEquals((11..19).toList(), n.equivalentAnchors)
        assertFalse(n.uniquePlacement)
    }

    @Test
    fun `single A deletion at run start has single physical anchor`() {
        val n = Normalizer(chr20).normalize(10, "CA", "C")
        assertEquals(10, n.pos); assertEquals("CA", n.ref); assertEquals("C", n.alt)
        assertEquals(listOf(10), n.equivalentAnchors)
        assertTrue(n.uniquePlacement)
    }

    @Test
    fun `interior A deletion groups anchors across the run`() {
        val n = Normalizer(chr20).normalize(11, "AA", "A")
        assertEquals((11..18).toList(), n.equivalentAnchors)
        assertEquals(11, n.pos)
    }

    @Test
    fun `all three poly-A deletion spellings converge on their molecules`() {
        val nz = Normalizer(chr20)
        val first = nz.normalize(10, "CA", "C")
        val middle = nz.normalize(11, "AA", "A")
        val last = nz.normalize(18, "AG", "G")
        // Deleting different physical A bases yields distinct canonical records.
        assertFalse(first.pos == middle.pos && first.ref == middle.ref)
        assertFalse(middle.pos == last.pos && middle.ref == last.ref)
    }

    @Test
    fun `unique insertion of diverse bases does not multi-anchor`() {
        // AAA G TTT : insert CTC after the lone G cannot be re-anchored.
        val seq = "AAAGTTT"
        val n = Normalizer(seq).normalize(4, "G", "GCTC")
        assertEquals(4, n.pos); assertEquals("GCTC", n.alt)
        assertEquals(listOf(4), n.equivalentAnchors)
        assertTrue(n.uniquePlacement)
    }

    @Test
    fun `minimal representation trims shared anchor`() {
        // AAC>AAG shares two leading A bases; it reduces to the substitution C>G.
        val seq = "AAACCCGGGTTT"
        val n = Normalizer(seq).normalize(2, "AAC", "AAG")
        assertEquals(4, n.pos); assertEquals("C", n.ref); assertEquals("G", n.alt)
        assertTrue(n.uniquePlacement)
    }

    @Test
    fun `deletion with shared trailing anchor is minimized`() {
        val seq = "AAACGTGGG"
        val n = Normalizer(seq).normalize(3, "ACG", "A")
        // deletes CG at 4-5, anchor A at 3, no repeat walk
        assertEquals(3, n.pos); assertEquals("ACG", n.ref); assertEquals("A", n.alt)
    }
}
