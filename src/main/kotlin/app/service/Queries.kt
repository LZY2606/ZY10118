package app.service

import app.model.CanonicalKey
import app.service.JsonParser.Companion.parse
import java.sql.ResultSet

/** Read-only queries used by the HTTP/UI layer. */
class Queries(private val svc: AppService) {

    fun listGroups(build: String?): List<Map<String, Any?>> {
        val sql = if (build == null)
            "SELECT * FROM canonical_groups ORDER BY build, contig, pos"
        else "SELECT * FROM canonical_groups WHERE build=? ORDER BY contig, pos"
        svc.db.connection.prepareStatement(sql).use { ps ->
            if (build != null) ps.setString(1, build)
            ps.executeQuery().use { rs ->
                val rows = ArrayList<Map<String, Any?>>()
                while (rs.next()) rows += groupRow(rs)
                return rows
            }
        }
    }

    fun groupDetail(canonicalId: String): Map<String, Any?> {
        svc.db.connection.prepareStatement("SELECT * FROM canonical_groups WHERE canonical_id=?").use {
            it.setString(1, canonicalId)
            it.executeQuery().use { rs ->
                if (!rs.next()) throw NotFoundException("unknown canonical group $canonicalId")
                val base = groupRow(rs)
                return base + mapOf(
                    "evidence" to evidence(canonicalId),
                    "projections" to projections(canonicalId),
                    "reviews" to reviews(canonicalId),
                )
            }
        }
    }

    fun evidence(canonicalId: String): List<Map<String, Any?>> {
        svc.db.connection.prepareStatement(
            """SELECT * FROM observations WHERE canonical_id=? OR
               (batch_id, line_number) IN (
                 SELECT batch_id, line_number FROM observations WHERE canonical_id=?
               ) ORDER BY batch_id, line_number, alt_index"""
        ).use {
            it.setString(1, canonicalId); it.setString(2, canonicalId)
            it.executeQuery().use { rs ->
                val out = ArrayList<Map<String, Any?>>()
                while (rs.next()) {
                    out += mapOf(
                        "batchId" to rs.getString("batch_id"),
                        "lineNumber" to rs.getInt("line_number"),
                        "contig" to rs.getString("contig"),
                        "pos" to rs.getInt("pos"),
                        "ref" to rs.getString("ref"),
                        "alt" to rs.getString("alt"),
                        "build" to rs.getString("build"),
                        "refStatus" to rs.getString("ref_status"),
                        "parentLineNumber" to rs.getObject("parent_line_number"),
                        "altIndex" to rs.getInt("alt_index"),
                        "childCount" to rs.getInt("child_count"),
                        "sourceId" to rs.getString("source_id"),
                        "canonicalId" to rs.getString("canonical_id"),
                    )
                }
                return out
            }
        }
    }

    fun projections(canonicalId: String): List<Map<String, Any?>> {
        svc.db.connection.prepareStatement(
            "SELECT * FROM projections WHERE canonical_id=? ORDER BY transcript_id"
        ).use {
            it.setString(1, canonicalId)
            it.executeQuery().use { rs ->
                val out = ArrayList<Map<String, Any?>>()
                while (rs.next()) {
                    out += mapOf(
                        "transcriptId" to rs.getString("transcript_id"),
                        "build" to rs.getString("build"),
                        "rulesVersion" to rs.getString("rules_version"),
                        "consequence" to rs.getString("consequence"),
                        "failure" to rs.getString("failure"),
                        "cdsStart" to rs.getObject("cds_start"),
                        "cdsEnd" to rs.getObject("cds_end"),
                        "proteinStart" to rs.getObject("protein_start"),
                        "proteinEnd" to rs.getObject("protein_end"),
                        "refCodon" to rs.getString("ref_codon"),
                        "altCodon" to rs.getString("alt_codon"),
                        "refAa" to rs.getString("ref_aa"),
                        "altAa" to rs.getString("alt_aa"),
                        "path" to (parse(rs.getString("path_json")) as List<*>),
                        "affectedExons" to (parse(rs.getString("affected_exons")) as List<*>),
                    )
                }
                return out
            }
        }
    }

    fun reviews(canonicalId: String): List<Map<String, Any?>> {
        svc.db.connection.prepareStatement(
            "SELECT * FROM reviews WHERE canonical_id=? ORDER BY version"
        ).use {
            it.setString(1, canonicalId)
            it.executeQuery().use { rs ->
                val out = ArrayList<Map<String, Any?>>()
                while (rs.next()) {
                    out += mapOf(
                        "version" to rs.getInt("version"),
                        "reviewer" to rs.getString("reviewer"),
                        "action" to rs.getString("action"),
                        "anchorOverride" to rs.getObject("anchor_override"),
                        "note" to rs.getString("note"),
                        "submittedAt" to rs.getString("submitted_at"),
                    )
                }
                return out
            }
        }
    }

    fun currentReviewVersion(canonicalId: String): Int =
        svc.db.connection.prepareStatement(
            "SELECT COALESCE(MAX(version),0) FROM reviews WHERE canonical_id=?"
        ).use {
            it.setString(1, canonicalId)
            it.executeQuery().use { r -> r.next(); r.getInt(1) }
        }

    fun listMismatches(build: String?): List<Map<String, Any?>> {
        val sql = "SELECT * FROM observations WHERE ref_status='MISMATCH'" +
            if (build == null) "" else " AND build=?"
        svc.db.connection.prepareStatement(sql).use { ps ->
            if (build != null) ps.setString(1, build)
            ps.executeQuery().use { rs ->
                val out = ArrayList<Map<String, Any?>>()
                while (rs.next()) out += mapOf(
                    "batchId" to rs.getString("batch_id"),
                    "lineNumber" to rs.getInt("line_number"),
                    "contig" to rs.getString("contig"),
                    "pos" to rs.getInt("pos"),
                    "ref" to rs.getString("ref"),
                    "alt" to rs.getString("alt"),
                    "build" to rs.getString("build"),
                    "sourceId" to rs.getString("source_id"),
                )
                return out
            }
        }
    }

    private fun groupRow(rs: ResultSet): Map<String, Any?> = mapOf(
        "canonicalId" to rs.getString("canonical_id"),
        "build" to rs.getString("build"),
        "contig" to rs.getString("contig"),
        "pos" to rs.getInt("pos"),
        "ref" to rs.getString("ref"),
        "alt" to rs.getString("alt"),
        "rulesVersion" to rs.getString("rules_version"),
        "equivalentAnchors" to (parse(rs.getString("equivalent_anchors")) as List<*>),
        "uniquePlacement" to (rs.getInt("unique_placement") == 1),
        "normSteps" to (parse(rs.getString("norm_steps")) as List<*>),
        "frozen" to (rs.getInt("frozen") == 1),
        "frozenAnchor" to rs.getObject("frozen_anchor"),
        "markedMismatch" to (rs.getInt("marked_mismatch") == 1),
    )

    fun batchReceipt(batchId: String): Map<String, Any?>? =
        svc.db.connection.prepareStatement("SELECT receipt FROM batches WHERE batch_id=?").use {
            it.setString(1, batchId)
            it.executeQuery().use { r ->
                if (!r.next()) null else (parse(r.getString(1)) as Map<String, Any?>)
            }
        }

    fun listBuilds(): List<Map<String, Any?>> {
        svc.db.connection.prepareStatement(
            "SELECT b.name, b.imported_at, (SELECT COUNT(*) FROM contigs c WHERE c.build=b.name) contigs FROM builds b"
        ).executeQuery().use { rs ->
            val out = ArrayList<Map<String, Any?>>()
            while (rs.next()) out += mapOf(
                "name" to rs.getString("name"),
                "importedAt" to rs.getString("imported_at"),
                "contigCount" to rs.getInt("contigs"),
            )
            return out
        }
    }

    fun alignmentWindow(canonicalId: String, radius: Int = 8): Map<String, Any?> {
        val g = groupDetail(canonicalId)
        val buildName = g["build"] as String
        val contig = g["contig"] as String
        val pos = g["pos"] as Int
        val seq = svc.build(buildName).sequence(contig)
            ?: throw NotFoundException("contig missing")
        val lo = maxOf(1, pos - radius)
        val hi = minOf(seq.length, pos + (g["ref"] as String).length - 1 + radius)
        return mapOf(
            "contig" to contig, "build" to buildName,
            "windowStart" to lo, "windowEnd" to hi,
            "referenceWindow" to seq.substring(lo - 1, hi),
            "alleleStart" to pos,
            "alleleEnd" to pos + (g["ref"] as String).length - 1,
        )
    }
}
