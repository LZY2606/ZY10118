package app.service

import app.model.*
import app.norm.Normalizer
import app.store.Database
import app.store.FastaReader
import app.transcript.Projector
import app.transcript.TranscriptLoader
import app.vcf.VcfParser
import java.io.File
import java.security.MessageDigest
import java.time.Instant

/** Central orchestration. All write paths are atomic via [Database.transaction]. */
class AppService(
    val db: Database,
    private val fixtureDir: File,
) {
    private val builds = HashMap<String, Build>()
    private val transcriptsByBuild = HashMap<String, List<Transcript>>()
    val queries = Queries(this)
    val reviews = Reviews(this)
    val mapImport = MapImport(this)

    fun loadBuild(name: String, fasta: File) {
        val build = FastaReader.read(fasta, name)
        synchronized(builds) { builds[name] = build }
        val now = Instant.now().toString()
        db.transaction {
            db.connection.prepareStatement(
                "INSERT OR IGNORE INTO builds(name, imported_at, fasta_path) VALUES(?,?,?)"
            ).use { it.setString(1, name); it.setString(2, now); it.setString(3, fasta.absolutePath); it.executeUpdate() }
            build.contigs.forEach { (contig, seq) ->
                db.connection.prepareStatement(
                    "INSERT OR IGNORE INTO contigs(build, contig, length, sequence) VALUES(?,?,?,?)"
                ).use {
                    it.setString(1, name); it.setString(2, contig)
                    it.setInt(3, seq.length); it.setString(4, seq); it.executeUpdate()
                }
            }
        }
    }

    fun loadTranscripts(name: String, tsv: File) {
        val rows = TranscriptLoader.read(tsv, name)
        synchronized(transcriptsByBuild) { transcriptsByBuild[name] = rows }
        db.transaction {
            rows.forEach { tx ->
                db.connection.prepareStatement(
                    """INSERT OR IGNORE INTO transcripts(build, transcript_id, contig, strand,
                      exon_starts, exon_ends, cds_start, cds_end, cds_exon_starts, cds_exon_ends)
                      VALUES(?,?,?,?,?,?,?,?,?,?)"""
                ).use { ps ->
                    ps.setString(1, tx.build); ps.setString(2, tx.id); ps.setString(3, tx.contig)
                    ps.setString(4, tx.strand)
                    ps.setString(5, tx.exonStarts.joinToString(","))
                    ps.setString(6, tx.exonEnds.joinToString(","))
                    ps.setInt(7, tx.cdsStart); ps.setInt(8, tx.cdsEnd)
                    ps.setString(9, tx.cdsExonStarts.joinToString(","))
                    ps.setString(10, tx.cdsExonEnds.joinToString(","))
                    ps.executeUpdate()
                }
            }
        }
    }

    fun build(name: String): Build =
        synchronized(builds) { builds[name] }
            ?: throw NotFoundException("build $name is not imported")

    fun transcripts(buildName: String): List<Transcript> =
        synchronized(transcriptsByBuild) { transcriptsByBuild[buildName] ?: emptyList() }

    /** Atomically import a VCF batch. Same batch ID retry is idempotent and
     *  returns the original receipt even if the payload differs. */
    fun importBatch(batchId: String, buildName: String, vcfContent: String): BatchReceipt {
        build(buildName) // fail early if unknown
        val sha = sha256(vcfContent)
        val existing = findReceipt(batchId)
        if (existing != null) {
            if (existing.getString("payload_sha256") != sha) {
                throw ConflictException(
                    "batch $batchId already exists with a different payload; the original receipt is retained"
                )
            }
            return receiptRow(existing)
        }
        val records = VcfParser.parse(vcfContent, batchId)
        return db.transaction {
            var accepted = 0
            var mismatches = 0
            val groupKeys = LinkedHashSet<String>()
            records.forEach { rec ->
                val seq = build(buildName).sequence(rec.contig)
                val refMatches = seq != null &&
                    rec.pos - 1 + rec.ref.length <= seq.length &&
                    seq.substring(rec.pos - 1, rec.pos - 1 + rec.ref.length)
                        .equals(rec.ref, ignoreCase = true)
                val multi = rec.alts.size > 1
                // Original composite evidence row (one row per source VCF line).
                insertObservation(rec, buildName, rec.ref, rec.alts.joinToString(","),
                    refMatches = refMatches, parent = null, altIndex = 0,
                    childCount = if (multi) rec.alts.size else 0)
                if (!refMatches) mismatches++
                rec.alts.forEachIndexed { idx, alt ->
                    // Split children preserve the parent provenance; they are
                    // separate alleles but never duplicate source evidence.
                    val childRefMatches = seq != null &&
                        rec.pos - 1 + rec.ref.length <= seq.length &&
                        seq.substring(rec.pos - 1, rec.pos - 1 + rec.ref.length)
                            .equals(rec.ref, ignoreCase = true)
                    insertObservation(rec, buildName, rec.ref, alt,
                        refMatches = childRefMatches,
                        parent = if (multi) rec.lineNumber else null,
                        altIndex = idx + 1,
                        childCount = 0)
                    if (childRefMatches) {
                        accepted++
                        val normalizer = Normalizer(seq!!)
                        val norm = normalizer.normalize(rec.pos, rec.ref, alt)
                        val key = CanonicalKey(buildName, rec.contig, norm.pos, norm.ref, norm.alt)
                        val cid = key.stableId()
                        groupKeys.add(cid)
                        upsertCanonical(buildName, norm, cid)
                        assignCanonical(rec.batchId, rec.lineNumber, idx + 1, cid)
                        projectAllele(buildName, key, norm, refMatches = true)
                    } else {
                        mismatches++
                    }
                }
            }
            val now = Instant.now().toString()
            val receipt = BatchReceipt(batchId, buildName, records.size, accepted,
                mismatches, groupKeys.size, duplicate = false, receivedAt = now)
            db.connection.prepareStatement(
                """INSERT INTO batches(batch_id, build, payload_sha256, received_at,
                   record_count, receipt) VALUES(?,?,?,?,?,?)"""
            ).use { ps ->
                ps.setString(1, batchId); ps.setString(2, buildName); ps.setString(3, sha)
                ps.setString(4, now); ps.setInt(5, records.size)
                ps.setString(6, Json.stringify(mapOf(
                    "batchId" to receipt.batchId, "build" to receipt.build,
                    "recordCount" to receipt.recordCount,
                    "acceptedAlleles" to receipt.acceptedAlleles,
                    "refMismatches" to receipt.refMismatches,
                    "groupsCreated" to receipt.groupsCreated,
                    "receivedAt" to receipt.receivedAt)))
                ps.executeUpdate()
            }
            receipt
        }
    }

    private fun insertObservation(
        rec: VcfRecord, build: String, ref: String, alt: String, refMatches: Boolean,
        parent: Int?, altIndex: Int, childCount: Int,
    ) {
        db.connection.prepareStatement(
            """INSERT INTO observations(batch_id, line_number, contig, pos, ref, alt,
               build, ref_status, parent_line_number, alt_index, child_count, source_id)
               VALUES(?,?,?,?,?,?,?,?,?,?,?,?)"""
        ).use { ps ->
            ps.setString(1, rec.batchId); ps.setInt(2, rec.lineNumber)
            ps.setString(3, rec.contig); ps.setInt(4, rec.pos)
            ps.setString(5, ref); ps.setString(6, alt); ps.setString(7, build)
            ps.setString(8, if (refMatches) "MATCH" else "MISMATCH")
            if (parent == null) ps.setNull(9, java.sql.Types.INTEGER) else ps.setInt(9, parent)
            ps.setInt(10, altIndex); ps.setInt(11, childCount)
            ps.setString(12, rec.ids.firstOrNull() ?: "line:${rec.lineNumber}")
            ps.executeUpdate()
        }
    }

    private fun assignCanonical(batchId: String, lineNumber: Int, altIndex: Int, cid: String) {
        db.connection.prepareStatement(
            "UPDATE observations SET canonical_id=? WHERE batch_id=? AND line_number=? AND alt_index=?"
        ).use {
            it.setString(1, cid); it.setString(2, batchId)
            it.setInt(3, lineNumber); it.setInt(4, altIndex); it.executeUpdate()
        }
    }

    private fun upsertCanonical(buildName: String, norm: Normalized, cid: String) {
        val steps = Json.stringify(norm.steps.map { step ->
            mapOf("op" to step.op, "pos" to step.pos, "ref" to step.ref,
                "alt" to step.alt, "note" to step.note)
        })
        db.connection.prepareStatement(
            """INSERT INTO canonical_groups(canonical_id, build, contig, pos, ref, alt,
               rules_version, equivalent_anchors, unique_placement, norm_steps)
               VALUES(?,?,?,?,?,?,?,?,?,?)
               ON CONFLICT(canonical_id) DO NOTHING"""
        ).use {
            it.setString(1, cid); it.setString(2, buildName); it.setString(3, "")
            // contig comes from the key but is not in Normalized; fill via update below
            it.setInt(4, norm.pos); it.setString(5, norm.ref); it.setString(6, norm.alt)
            it.setString(7, ANNOTATION_RULES_VERSION)
            it.setString(8, Json.stringify(norm.equivalentAnchors))
            it.setInt(9, if (norm.uniquePlacement) 1 else 0)
            it.setString(10, steps); it.executeUpdate()
        }
        val contig = cid.substringAfter("|").substringBefore("|")
        db.connection.prepareStatement("UPDATE canonical_groups SET contig=? WHERE canonical_id=? AND contig=''")
            .use { it.setString(1, contig); it.setString(2, cid); it.executeUpdate() }
    }

    private fun projectAllele(buildName: String, key: CanonicalKey, norm: Normalized,
                              refMatches: Boolean) {
        val projector = Projector(buildName, build(buildName).contigs)
        transcripts(buildName).filter { it.contig == key.contig }.forEach { tx ->
            val p = projector.project(tx, norm, refMatches)
            saveProjection(key.stableId(), p)
        }
    }

    private fun saveProjection(cid: String, p: Projection) {
        db.connection.prepareStatement(
            """INSERT INTO projections(canonical_id, transcript_id, build, rules_version,
               consequence, failure, cds_start, cds_end, protein_start, protein_end,
               ref_codon, alt_codon, ref_aa, alt_aa, path_json, affected_exons)
               VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
               ON CONFLICT(canonical_id, transcript_id) DO NOTHING"""
        ).use { ps ->
            ps.setString(1, cid); ps.setString(2, p.transcriptId); ps.setString(3, p.build)
            ps.setString(4, p.rulesVersion); ps.setString(5, p.consequence)
            ps.setString(6, p.failure)
            setIntOrNull(ps, 7, p.cdsStart); setIntOrNull(ps, 8, p.cdsEnd)
            setIntOrNull(ps, 9, p.proteinStart); setIntOrNull(ps, 10, p.proteinEnd)
            ps.setString(11, p.refCodon); ps.setString(12, p.altCodon)
            ps.setString(13, p.refAa); ps.setString(14, p.altAa)
            ps.setString(15, Json.stringify(p.path))
            ps.setString(16, Json.stringify(p.affectedExons))
            ps.executeUpdate()
        }
    }

    private fun setIntOrNull(ps: java.sql.PreparedStatement, idx: Int, v: Int?) {
        if (v == null) ps.setNull(idx, java.sql.Types.INTEGER) else ps.setInt(idx, v)
    }

    private fun findReceipt(batchId: String): java.sql.ResultSet? {
        val ps = db.connection.prepareStatement("SELECT * FROM batches WHERE batch_id=?")
        ps.setString(1, batchId)
        val rs = ps.executeQuery()
        return if (rs.next()) rs else { rs.close(); ps.close(); null }
    }

    private fun receiptRow(rs: java.sql.ResultSet): BatchReceipt = BatchReceipt(
        rs.getString("batch_id"), rs.getString("build"), rs.getInt("record_count"),
        acceptedAlleles = countObs(rs.getString("batch_id"), "MATCH", child = true),
        refMismatches = countObs(rs.getString("batch_id"), "MISMATCH", child = false),
        groupsCreated = countGroupsForBatch(rs.getString("batch_id")),
        duplicate = true,
        receivedAt = rs.getString("received_at"),
    )

    private fun countObs(batchId: String, status: String, child: Boolean): Int =
        db.connection.prepareStatement(
            "SELECT COUNT(*) FROM observations WHERE batch_id=? AND ref_status=? AND alt_index>0"
        ).use {
            it.setString(1, batchId); it.setString(2, status)
            it.executeQuery().use { r -> r.next(); r.getInt(1) }
        }

    private fun countGroupsForBatch(batchId: String): Int = db.connection.prepareStatement(
        "SELECT COUNT(DISTINCT canonical_id) FROM observations WHERE batch_id=? AND canonical_id IS NOT NULL"
    ).use {
        it.setString(1, batchId); it.executeQuery().use { r -> r.next(); r.getInt(1) }
    }

    private fun sha256(text: String): String =
        MessageDigest.getInstance("SHA-256").digest(text.toByteArray())
            .joinToString("") { "%02x".format(it) }

    private fun nowIso(): String = Instant.now().toString()
}
