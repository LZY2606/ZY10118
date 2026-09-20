package app

import app.norm.Normalizer
import app.transcript.Projector
import app.transcript.FailureType
import org.junit.jupiter.api.Assertions.assertEquals
import org.junit.jupiter.api.Assertions.assertNull
import org.junit.jupiter.api.Assertions.assertTrue
import org.junit.jupiter.api.Test

class ProjectorTest {
    private val build = FixtureSupport.b37()
    private val txBeta = FixtureSupport.transcriptsB37().first { it.id == "TX_BETA" }
    private val txDelta = FixtureSupport.transcriptsB37().first { it.id == "TX_DELTA" }
    private val projector = Projector("b37", build.contigs)

    private fun norm(pos: Int, ref: String, alt: String) =
        Normalizer(build.sequence("chr11")!!).normalize(pos, ref, alt)

    @Test
    fun `start codon substitution is START_LOST on plus strand`() {
        val allele = norm(40, "A", "C")
        val p = projector.project(txBeta, allele, refMatches = true)
        assertEquals("START_LOST", p.consequence)
        assertNull(p.failure)
        assertEquals(1, p.cdsStart)
        assertEquals("ATG", p.refCodon); assertEquals("CTG", p.altCodon)
        assertEquals("M", p.refAa); assertEquals("L", p.altAa)
        assertEquals(1, p.proteinStart)
    }

    @Test
    fun `stop codon disruption is STOP_LOST`() {
        val allele = norm(167, "T", "A")
        val p = projector.project(txBeta, allele, true)
        assertEquals("STOP_LOST", p.consequence)
        assertEquals("TAA", p.refCodon); assertEquals("AAA", p.altCodon)
    }

    @Test
    fun `codon straddling exon boundary uses spliced bases`() {
        // g79,g80 (exon1) + g121 (exon2) form codon 14 = TTA.
        val allele = norm(80, "T", "C")
        val p = projector.project(txBeta, allele, true)
        assertEquals("MISSENSE_VARIANT", p.consequence)
        assertEquals(14, p.proteinStart)
        assertEquals("TTA", p.refCodon); assertEquals("TCA", p.altCodon)
        assertTrue(p.path.any { it.contains("codon 14") })
    }

    @Test
    fun `intronic base is INTRON_VARIANT`() {
        val allele = norm(100, "A", "T")
        val p = projector.project(txBeta, allele, true)
        assertEquals("INTRON_VARIANT", p.consequence)
        assertNull(p.failure)
    }

    @Test
    fun `five prime UTR base is classified`() {
        val allele = norm(30, "C", "G")
        val p = projector.project(txBeta, allele, true)
        assertEquals("FIVE_PRIME_UTR_VARIANT", p.consequence)
    }

    @Test
    fun `far intergenic allele fails OUTSIDE_TRANSCRIPT`() {
        val allele = norm(5, "G", "T")
        val p = projector.project(txBeta, allele, true)
        assertEquals("PROJECTION_FAILED", p.consequence)
        assertEquals(FailureType.OUTSIDE_TRANSCRIPT, p.failure)
    }

    @Test
    fun `unique in-frame deletion is INFRAME_INDEL`() {
        val allele = Normalizer(build.sequence("chr11")!!).normalize(160, "GACT", "G")
        assertEquals(listOf(160), allele.equivalentAnchors)
        val p = projector.project(txBeta, allele, true)
        assertEquals("INFRAME_INDEL", p.consequence)
    }

    @Test
    fun `minus strand start codon maps correctly`() {
        // g280 is the first coding base read on the minus strand.
        val allele = Normalizer(build.sequence("chr17")!!).normalize(280, "T", "C")
        val p = Projector("b37", build.contigs).project(txDelta, allele, true)
        assertEquals("START_LOST", p.consequence)
        assertEquals(1, p.cdsStart)
        assertEquals("ATG", p.refCodon); assertEquals("CTG", p.altCodon)
        assertEquals(2, (p.affectedExons.firstOrNull() ?: -1))
    }

    @Test
    fun `multi-anchor allele withholds projection`() {
        val allele = Normalizer(build.sequence("chr20")!!).normalize(14, "A", "AA")
        // No chr20 transcript in fixtures; construct a throwaway expectation via
        // the BETA transcript contig mismatch is not relevant, so test the gate
        // directly on chr20 with a projector that has the sequence.
        val probe = app.model.Transcript(
            "PROBE", "b37", "chr20", "+", listOf(1), listOf(130),
            1, 130, listOf(1), listOf(130),
        )
        val p = Projector("b37", build.contigs).project(probe, allele, true)
        assertEquals(FailureType.NON_UNIQUE_LEFT_ALIGN, p.failure)
    }
}
