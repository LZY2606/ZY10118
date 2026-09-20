package app.service

import app.model.CanonicalKey
import app.norm.Normalizer
import java.io.File
import java.time.Instant

/** Review actions with optimistic concurrency. */
class Reviews(private val svc: AppService) {

    /** Submit a review. [expectedVersion] must equal the current stored version;
     *  otherwise a 409 conflict is returned and nothing is overwritten. */
    fun submit(
        canonicalId: String, reviewer: String, action: String,
        expectedVersion: Int, anchorOverride: Int?, note: String?,
    ): Map<String, Any?> {
        require(reviewer.isNotBlank()) { "reviewer is required" }
        require(action in setOf("FREEZE", "MARK_MISMATCH", "CLEAR", "COMMENT")) {
            "unknown action $action"
        }
        return svc.db.transaction {
            val current = svc.queries.currentReviewVersion(canonicalId)
            if (current != expectedVersion) {
                throw ConflictException(
                    "canonical $canonicalId is at version $current; expected $expectedVersion"
                )
            }
            val newVersion = current + 1
            val now = Instant.now().toString()
            svc.db.connection.prepareStatement(
                """INSERT INTO reviews(canonical_id, version, reviewer, action,
                   anchor_override, note, submitted_at) VALUES(?,?,?,?,?,?,?)"""
            ).use { ps ->
                ps.setString(1, canonicalId); ps.setInt(2, newVersion)
                ps.setString(3, reviewer); ps.setString(4, action)
                if (anchorOverride == null) ps.setNull(5, java.sql.Types.INTEGER)
                else ps.setInt(5, anchorOverride)
                ps.setString(6, note); ps.setString(7, now)
                ps.executeUpdate()
            }
            when (action) {
                "FREEZE" -> {
                    val anchor = anchorOverride ?: defaultAnchor(canonicalId)
                    svc.db.connection.prepareStatement(
                        "UPDATE canonical_groups SET frozen=1, frozen_anchor=?, marked_mismatch=0 WHERE canonical_id=?"
                    ).use { it.setInt(1, anchor); it.setString(2, canonicalId); it.executeUpdate() }
                }
                "MARK_MISMATCH" -> svc.db.connection.prepareStatement(
                    "UPDATE canonical_groups SET marked_mismatch=1, frozen=0, frozen_anchor=NULL WHERE canonical_id=?"
                ).use { it.setString(1, canonicalId); it.executeUpdate() }
                "CLEAR" -> svc.db.connection.prepareStatement(
                    "UPDATE canonical_groups SET frozen=0, frozen_anchor=NULL, marked_mismatch=0 WHERE canonical_id=?"
                ).use { it.setString(1, canonicalId); it.executeUpdate() }
            }
            mapOf("canonicalId" to canonicalId, "version" to newVersion,
                "reviewer" to reviewer, "action" to action, "submittedAt" to now)
        }
    }

    private fun defaultAnchor(canonicalId: String): Int {
        svc.db.connection.prepareStatement(
            "SELECT pos FROM canonical_groups WHERE canonical_id=?"
        ).use {
            it.setString(1, canonicalId)
            it.executeQuery().use { rs ->
                if (!rs.next()) throw NotFoundException("unknown group $canonicalId")
                return rs.getInt("pos")
            }
        }
    }
}

/** Import a build-to-build chain and produce an impact report. Atomic: the
 *  segments and every impact row commit together or not at all. */
class MapImport(private val svc: AppService) {
    fun importChain(
        sourceBuild: String, destBuild: String, chainTsv: String,
    ): ImpactReport {
        val dest = svc.build(destBuild)
        val src = svc.build(sourceBuild)
        val sha = java.security.MessageDigest.getInstance("SHA-256")
            .digest(chainTsv.toByteArray()).joinToString("") { "%02x".format(it) }
        val segments = app.liftover.LiftOver.let {
            val tmp = File.createTempFile("chain", ".tsv")
            tmp.deleteOnExit(); tmp.writeText(chainTsv)
            app.liftover.LiftOver.parse(tmp, sourceBuild, destBuild).segments
        }
        return svc.db.transaction {
            val liftoverId: Int
            svc.db.connection.prepareStatement(
                """INSERT INTO liftovers(source_build, dest_build, imported_at, payload_sha256)
                   VALUES(?,?,?,?)
                   ON CONFLICT(source_build, dest_build) DO UPDATE SET imported_at=excluded.imported_at,
                     payload_sha256=excluded.payload_sha256"""
            ).use { ps ->
                ps.setString(1, sourceBuild); ps.setString(2, destBuild)
                ps.setString(3, Instant.now().toString()); ps.setString(4, sha)
                ps.executeUpdate()
            }
            svc.db.connection.prepareStatement(
                "DELETE FROM liftover_segments WHERE liftover_id=(SELECT id FROM liftovers WHERE source_build=? AND dest_build=?)"
            ).use { it.setString(1, sourceBuild); it.setString(2, destBuild); it.executeUpdate() }
            svc.db.connection.prepareStatement(
                "SELECT id FROM liftovers WHERE source_build=? AND dest_build=?"
            ).use {
                it.setString(1, sourceBuild); it.setString(2, destBuild)
                it.executeQuery().use { r -> r.next(); liftoverId = r.getInt(1) }
            }
            segments.forEachIndexed { idx, seg ->
                svc.db.connection.prepareStatement(
                    """INSERT INTO liftover_segments(liftover_id, ord, src_contig, src_start,
                       src_end, dst_contig, dst_start, dst_end, dst_strand)
                       VALUES(?,?,?,?,?,?,?,?,?)"""
                ).use { ps ->
                    ps.setInt(1, liftoverId); ps.setInt(2, idx + 1)
                    ps.setString(3, seg.srcContig); ps.setInt(4, seg.srcStart)
                    ps.setInt(5, seg.srcEnd); ps.setString(6, seg.dstContig)
                    ps.setInt(7, seg.dstStart); ps.setInt(8, seg.dstEnd)
                    ps.setString(9, seg.dstStrand); ps.executeUpdate()
                }
            }

            val chain = app.liftover.LiftOver(sourceBuild, destBuild, segments)
            val items = ArrayList<ImpactItem>()
            svc.db.connection.prepareStatement(
                "SELECT * FROM canonical_groups WHERE build=?"
            ).use {
                it.setString(1, sourceBuild)
                it.executeQuery().use { rs ->
                    while (rs.next()) {
                        val norm = app.model.Normalized(
                            pos = rs.getInt("pos"),
                            ref = rs.getString("ref"), alt = rs.getString("alt"),
                            steps = emptyList(),
                            equivalentAnchors = (JsonParser.parse(
                                rs.getString("equivalent_anchors")) as List<*>).map { v -> (v as Number).toInt() },
                            uniquePlacement = rs.getInt("unique_placement") == 1,
                        )
                        val cid = rs.getString("canonical_id")
                        val contig = rs.getString("contig")
                        val lr = chain.lift(contig, norm, dest.contigs)
                        items += ImpactItem(cid, lr.status, lr.destContig, lr.destPos,
                            lr.segmentsUsed, lr.note)
                    }
                }
            }
            svc.db.connection.prepareStatement(
                "DELETE FROM impact_reports WHERE liftover_id=?"
            ).use { it.setInt(1, liftoverId); it.executeUpdate() }
            items.forEach { item ->
                svc.db.connection.prepareStatement(
                    """INSERT INTO impact_reports(liftover_id, canonical_id, status,
                       dest_contig, dest_pos, dest_ref, dest_alt, segments_used, note, created_at)
                       VALUES(?,?,?,?,?,?,?,?,?,?)"""
                ).use { ps ->
                    ps.setInt(1, liftoverId); ps.setString(2, item.canonicalId)
                    ps.setString(3, item.status)
                    ps.setString(4, item.destContig)
                    if (item.destPos == null) ps.setNull(5, java.sql.Types.INTEGER)
                    else ps.setInt(5, item.destPos)
                    ps.setNull(6, java.sql.Types.VARCHAR); ps.setNull(7, java.sql.Types.VARCHAR)
                    ps.setString(8, Json.stringify(item.segments))
                    ps.setString(9, item.note); ps.setString(10, Instant.now().toString())
                    ps.executeUpdate()
                }
            }
            val summary = items.groupingBy { it.status }.eachCount()
            ImpactReport(sourceBuild, destBuild, items.size, items, summary)
        }
    }
}
