package app

import app.service.ConflictException
import app.service.AppService
import app.store.Database
import app.transcript.FailureType
import org.junit.jupiter.api.AfterEach
import org.junit.jupiter.api.Assertions.assertEquals
import org.junit.jupiter.api.Assertions.assertNotEquals
import org.junit.jupiter.api.Assertions.assertNotNull
import org.junit.jupiter.api.Assertions.assertThrows
import org.junit.jupiter.api.Assertions.assertTrue
import org.junit.jupiter.api.BeforeEach
import org.junit.jupiter.api.Test
import org.junit.jupiter.api.io.TempDir
import java.io.File

class ServiceIntegrationTest {
    @TempDir lateinit var tmp: File
    private lateinit var db: Database
    private lateinit var svc: AppService

    @BeforeEach
    fun setup() {
        db = Database(File(tmp, "state.db").absolutePath)
        svc = AppService(db, FixtureSupport.dir)
        svc.loadBuild("b37", File(FixtureSupport.dir, "ref_b37.fa"))
        svc.loadBuild("b38", File(FixtureSupport.dir, "ref_b38.fa"))
        svc.loadTranscripts("b37", File(FixtureSupport.dir, "transcripts_b37.tsv"))
        svc.loadTranscripts("b38", File(FixtureSupport.dir, "transcripts_b38.tsv"))
    }

    @AfterEach
    fun tearDown() { db.close() }

    @Test
    fun `seed import creates expected groups and projections`() {
        val receipt = svc.importBatch("B1", "b37", FixtureSupport.seedVcf())
        assertEquals(15, receipt.recordCount)
        val groups = svc.queries.listGroups("b37")
        // 15 VCF rows, one multi-allelic row contributes two children -> groups:
        // at least the documented unique plus poly-A groups.
        assertTrue(groups.isNotEmpty())
        // start-codon group carries a START_LOST projection.
        val start = svc.queries.groupDetail("b37|chr11|40|A|C")
        val beta = start["projections"] as List<*>
        assertTrue(beta.any { (it as Map<*, *>)["consequence"] == "START_LOST" })
    }

    @Test
    fun `identical batch id retry returns original receipt`() {
        val first = svc.importBatch("DUP", "b37", FixtureSupport.seedVcf())
        val second = svc.importBatch("DUP", "b37", FixtureSupport.seedVcf())
        assertEquals(first.receivedAt, second.receivedAt)
        assertEquals(true, second.duplicate)
        // importing the same id with a different payload is rejected
        assertThrows(ConflictException::class.java) {
            svc.importBatch("DUP", "b37", "chr20\t1\t.\tA\tT\n")
        }
    }

    @Test
    fun `split children trace parent without duplicating source evidence`() {
        svc.importBatch("MULTI", "b37",
            "chr20\t14\tX\tA\tAA,C\n")
        val groups = svc.queries.listGroups("b37").associateBy { it["canonicalId"] }
        // both children normalized under distinct keys (insertion vs substitution)
        assertNotNull(groups["b37|chr20|11|A|AA"])
        assertNotNull(groups["b37|chr20|14|A|C"])
        val insertion = svc.queries.groupDetail("b37|chr20|11|A|AA")
        val evidence = insertion["evidence"] as List<Map<String, Any?>>
        // exactly one parent row and at least one child for the same source line
        val parents = evidence.count { it["altIndex"] == 0 }
        assertEquals(1, parents)
        assertTrue(evidence.any { it["parentLineNumber"] != null })
    }

    @Test
    fun `reference mismatch is recorded and projection withheld`() {
        svc.importBatch("MM", "b37", "chr20\t25\tBAD\tZZZ\tZ\n")
        val mismatches = svc.queries.listMismatches("b37")
        assertTrue(mismatches.any { it["sourceId"] == "BAD" })
        // no canonical group should have been created for the mismatched allele
        assertTrue(svc.queries.listGroups("b37").none { it["pos"] == 25 })
    }

    @Test
    fun `optimistic review lock blocks stale submissions`() {
        svc.importBatch("R", "b37", FixtureSupport.seedVcf())
        val cid = "b37|chr11|40|A|C"
        svc.reviews.submit(cid, "alice", "FREEZE", 0, null, null)
        assertThrows(ConflictException::class.java) {
            svc.reviews.submit(cid, "mallory", "MARK_MISMATCH", 0, null, "too late")
        }
        // newer version proceeds
        svc.reviews.submit(cid, "bob", "MARK_MISMATCH", 1, null, "ok")
        val detail = svc.queries.groupDetail(cid)
        assertEquals(true, detail["markedMismatch"])
        assertEquals(listOf(1, 2), (detail["reviews"] as List<Map<String, Any?>>).map { it["version"] })
    }

    @Test
    fun `chain import reports all liftover states`() {
        svc.importBatch("L", "b37", FixtureSupport.seedVcf())
        val report = svc.mapImport.importChain("b37", "b38", FixtureSupport.chain())
        val statuses = report.items.map { it.status }.toSet()
        assertTrue(statuses.contains("MAPPED"))
        assertTrue(statuses.contains("COLLISION"))
        assertTrue(statuses.contains("STRAND_FLIPPED"))
        assertTrue(statuses.contains("GAP"))
    }

    @Test
    fun `old b37 projections remain after b38 chain import`() {
        svc.importBatch("OLD", "b37", FixtureSupport.seedVcf())
        svc.mapImport.importChain("b37", "b38", FixtureSupport.chain())
        val detail = svc.queries.groupDetail("b37|chr11|40|A|C")
        val projections = detail["projections"] as List<Map<String, Any?>>
        // b37 result is preserved and still tagged with build b37
        assertTrue(projections.any { it["build"] == "b37" && it["consequence"] == "START_LOST" })
    }

    @Test
    fun `failed batch rolls back atomically`() {
        // first record valid, second structurally invalid -> whole batch rejected
        val bad = "chr20\t1\t.\tA\tT\nchr20\tX\t.\tA\tT\n"
        assertThrows(Exception::class.java) {
            svc.importBatch("ATOMIC", "b37", bad)
        }
        assertTrue(svc.queries.batchReceipt("ATOMIC") == null)
        assertTrue(svc.queries.listGroups("b37").none { true && it["pos"] == 1 })
    }
}
