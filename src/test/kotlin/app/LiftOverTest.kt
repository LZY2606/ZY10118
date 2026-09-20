package app

import app.liftover.LiftOver
import app.norm.Normalizer
import org.junit.jupiter.api.Assertions.assertEquals
import org.junit.jupiter.api.Assertions.assertTrue
import org.junit.jupiter.api.Test
import java.io.File

class LiftOverTest {
    private val b37 = FixtureSupport.b37()
    private val b38 = FixtureSupport.b38()
    private val chain = LiftOver.parse(File(FixtureSupport.dir, "map_b37_to_b38.tsv"), "b37", "b38")

    @Test
    fun `clean offset chain maps with shift`() {
        val allele = Normalizer(b37.sequence("chr20")!!).normalize(1, "A", "T")
        val r = chain.lift("chr20", allele, b38.contigs)
        assertEquals("MAPPED", r.status)
        assertEquals(6, r.destPos)   // chr20 1-50 -> b38 6-55
        assertEquals(3, r.segmentsUsed.single())
    }

    @Test
    fun `reverse strand chain flips and reverse complements`() {
        val allele = Normalizer(b37.sequence("chr17")!!).normalize(280, "T", "C")
        val r = chain.lift("chr17", allele, b38.contigs)
        assertEquals("STRAND_FLIPPED", r.status)
        assertEquals(21, r.destPos)   // 300 - 280 + 1
        // reverse complement of REF T -> A, ALT C -> G
        assertEquals("A", r.destRef); assertEquals("G", r.destAlt)
    }

    @Test
    fun `allele spanning adjacent segments is GAP not collision`() {
        // g48 CGAC>C spans src 50|51 break.
        val allele = Normalizer(b37.sequence("chr20")!!).normalize(48, "CGAC", "C")
        val r = chain.lift("chr20", allele, b38.contigs)
        assertEquals("GAP", r.status)
        assertTrue(r.segmentsUsed.size >= 2)
    }

    @Test
    fun `overlapping second chain is many-to-one COLLISION`() {
        val allele = Normalizer(b37.sequence("chr20")!!).normalize(20, "A", "T")
        val r = chain.lift("chr20", allele, b38.contigs)
        assertEquals("COLLISION", r.status)
        assertTrue(r.segmentsUsed.size >= 2)
    }

    @Test
    fun `plus strand chr11 applies plus ten offset`() {
        val allele = Normalizer(b37.sequence("chr11")!!).normalize(40, "A", "C")
        val r = chain.lift("chr11", allele, b38.contigs)
        assertEquals("MAPPED", r.status)
        assertEquals(50, r.destPos)
    }
}
