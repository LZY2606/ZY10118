package app

import org.junit.jupiter.api.Assertions.*
import org.junit.jupiter.api.Test

class IntegrationTest {
    private fun harness() = FixtureHarness.newService(FixtureHarness.tempDir())

    private fun vcf(vararg body: String) =
        listOf("##fileformat=VCFv4.2",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
            .plus(body).joinToString("\n") + "\n"
    private fun rec(chrom: String, pos: Int, ref: String, alt: String, src: String) =
        "$chrom\t$pos\t.\t$ref\t$alt\t.\tPASS\tSRC=$src"

    @Test
    fun `batch import is atomic and retry returns the original receipt`() {
        val (svc, vcf, _) = harness()
        val first = svc.importBatch("B1", "hg19", vcf)
        assertFalse(first.replay)
        assertEquals(10, first.records)
        assertEquals(11, first.alleles) // one composite line carries 2 ALTs
        assertEquals(1, first.refMismatches)

        val second = svc.importBatch("B1", "hg19", vcf)
        assertTrue(second.replay)
        assertEquals(first.receivedAt, second.receivedAt)
        assertEquals(first.groups, second.groups)

        // replay must not create duplicate raw evidence
        val groups = svc.listGroups("hg19")
        assertTrue(groups.isNotEmpty())
        groups.forEach { g ->
            val detail = svc.getGroup(g.id)!!
            val rawIds = detail.evidence.map { it.rawRecordId }.distinct()
            val compositeOriginal = detail.evidence.count { it.role == "ORIGINAL" }
            // each underlying raw record contributes at most one ORIGINAL row
            detail.evidence.groupBy { it.rawRecordId }.values.forEach { rows ->
                assertTrue(rows.count { it.role == "ORIGINAL" } <= 1)
            }
            assertTrue(rawIds.size >= compositeOriginal)
        }
    }

    @Test
    fun `replaying same batch id against different build is rejected`() {
        val (svc, vcf, _) = harness()
        svc.importBatch("B9", "hg19", vcf)
        val ex = assertThrows(IllegalArgumentException::class.java) {
            svc.importBatch("B9", "hg38", vcf)
        }
        assertTrue(ex.message!!.contains("already imported"))
    }

    @Test
    fun `composite line and split alleles are bidirectionally traceable`() {
        val (svc, vcf, _) = harness()
        svc.importBatch("B2", "hg19", vcf)
        val groups = svc.listGroups("hg19").map { svc.getGroup(it.id)!! }
        val compositeGroups = groups.filter { g ->
            g.evidence.any { it.altField.contains(",") }
        }
        assertEquals(2, compositeGroups.size) // two split alleles from the one composite line
        compositeGroups.forEach { g ->
            val rawIds = g.evidence.filter { it.altField.contains(",") }.map { it.rawRecordId }.distinct()
            rawIds.forEach { rid ->
                assertTrue(g.evidence.any { it.rawRecordId == rid && it.role == "DERIVED" })
                val originals = g.evidence.filter { it.rawRecordId == rid && it.role == "ORIGINAL" }
                assertEquals(1, originals.size) { "split replay must not duplicate original evidence for raw $rid" }
            }
        }
    }

    @Test
    fun `same biological change from different batches merges and keeps both provenance`() {
        val (svc, _, _) = harness()
        // Two representations of the same CAT deletion inside the CATCATCAT run at 899..909
        svc.importBatch("P1", "hg19", vcf(rec("chr1", 900, "CATC", "C", "labX")))
        svc.importBatch("P2", "hg19", vcf(rec("chr1", 903, "CATC", "C", "labY")))
        val matching = svc.listGroups("hg19").map { svc.getGroup(it.id)!! }
            .filter { it.evidence.any { e -> e.batchId == "P1" } && it.evidence.any { e -> e.batchId == "P2" } }
        assertEquals(1, matching.size)
        val g = matching.single()
        val batches = g.evidence.map { it.batchId }.distinct().toSet()
        assertEquals(setOf("P1", "P2"), batches)
    }

    @Test
    fun `projections cover exon, intron, cross-exon, negative strand and out-of-transcript`() {
        val (svc, vcf, _) = harness()
        svc.importBatch("B3", "hg19", vcf)
        val all = svc.listGroups("hg19").map { svc.getGroup(it.id)!! }

        // exonic SNV on positive transcript
        val exonic = all.first { g -> g.consequences.any { it.region == "EXONIC" && it.txId == "NM_DEMO_P" } }
        val cP = exonic.consequences.first { it.txId == "NM_DEMO_P" }
        assertEquals("PROJECTED", cP.status)
        assertNotNull(cP.hgvsgLike)
        assertTrue(cP.pathJson.contains("positive strand"))

        // negative-strand transcript gets reverse complemented allele text
        val neg = all.flatMap { it.consequences }.first { it.txId == "NM_DEMO_N" && it.status == "PROJECTED" && it.region == "EXONIC" }
        assertTrue(neg.pathJson.contains("reverse-complemented"))

        // cross-exon-boundary deletion -> EXON_STRADDLE failure on NM_DEMO_P
        val straddle = all.flatMap { it.consequences }
            .first { it.txId == "NM_DEMO_P" && it.failureType == "EXON_STRADDLE" }
        assertTrue(straddle.pathJson.contains("intron"))

        // a group far from both transcripts -> OUT_OF_TRANSCRIPT
        assertTrue(all.flatMap { it.consequences }.any {
            it.failureType == "OUT_OF_TRANSCRIPT"
        })
    }

    @Test
    fun `review uses optimistic locking and rejects stale versions`() {
        val (svc, vcf, _) = harness()
        svc.importBatch("B4", "hg19", vcf)
        val gid = svc.listGroups("hg19").first().id
        val v1 = svc.submitReview(gid, 0, "alice", "FROZEN", "ok")
        assertEquals(1L, v1)
        val v2 = svc.submitReview(gid, 1, "bob", "REF_FLAGGED", "new look")
        assertEquals(2L, v2)
        val ex = assertThrows(IllegalArgumentException::class.java) {
            svc.submitReview(gid, 1, "carol", "FROZEN", "stale")
        }
        assertTrue(ex.message!!.contains("version 2"))
        assertEquals("REF_FLAGGED", svc.getGroup(gid)!!.reviewState)
    }

    @Test
    fun `liftover import reports impact across unique flipped many-to-one interrupted unmapped`() {
        val (svc, vcf, chain) = harness()
        svc.importBatch("B5", "hg19", vcf)
        val impact = svc.importChain(chain)
        assertEquals("hg19", impact.fromBuild)
        assertEquals("hg38", impact.toBuild)
        assertTrue(impact.newUniqueGroups >= 1)
        assertTrue(impact.flippedGroups + impact.manyToOneGroups +
            impact.interruptedGroups + impact.unmappedGroups >= 1)

        val statuses = svc.listGroups("hg19").map { svc.getGroup(it.id)!! }
            .flatMap { it.liftovers }.map { it.status }.toSet()
        // the fixture chain exercises every named state among existing groups
        listOf("UNIQUE", "FLIPPED", "MANY_TO_ONE", "INTERRUPTED", "UNMAPPED").forEach {
            assertTrue(statuses.contains(it), "expected liftover status $it, got $statuses")
        }
    }

    @Test
    fun `old build results stay reproducible after a chain import`() {
        val (svc, vcf, chain) = harness()
        svc.importBatch("B6", "hg19", vcf)
        val before = svc.listGroups("hg19").map { svc.getGroup(it.id)!! }.map {
            it.id to (it.canonPos to it.canonRef)
        }.toMap()
        svc.importChain(chain)
        val after = svc.listGroups("hg19").map { svc.getGroup(it.id)!! }.map {
            it.id to (it.canonPos to it.canonRef)
        }.toMap()
        assertEquals(before, after)
    }

    @Test
    fun `ref mismatch group is flagged and carries expected reference`() {
        val (svc, vcf, _) = harness()
        svc.importBatch("B7", "hg19", vcf)
        val bad = svc.listGroups("hg19").map { svc.getGroup(it.id)!! }.first { it.refMismatch }
        assertNotNull(bad.expectedRef)
        assertTrue(bad.normSteps.any { it.kind == "ref_check" })
        assertTrue(bad.placements.isEmpty())
        assertTrue(bad.consequences.none { it.status == "PROJECTED" && it.region != null })
    }
}
