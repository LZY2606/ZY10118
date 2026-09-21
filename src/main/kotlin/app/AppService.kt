package app

import java.nio.file.Path
import java.security.MessageDigest
import java.time.Instant
import java.sql.ResultSet

data class GroupDetail(
    val id: Long,
    val build: String,
    val chrom: String,
    val canonPos: Int?,
    val canonRef: String?,
    val canonAlt: String?,
    val refMismatch: Boolean,
    val expectedRef: String?,
    val placements: List<Placement>,
    val normSteps: List<NormStep>,
    val preferenceRule: String,
    val rulesVersion: String,
    val frozen: Boolean,
    val reviewState: String,
    val evidence: List<EvidenceRow>,
    val consequences: List<ConsequenceRow>,
    val liftovers: List<LiftoverRowDto>,
    val reviewVersion: Long,
)

data class EvidenceRow(
    val rawRecordId: Long,
    val batchId: String,
    val source: String,
    val lineNumber: Int,
    val chrom: String,
    val pos: Int,
    val ref: String,
    val altField: String,
    val altIndex: Int,
    val role: String,
    val rawLine: String,
    val info: String,
)

data class LiftoverRowDto(
    val toBuild: String,
    val status: String,
    val dstChrom: String?,
    val dstPos: Int?,
    val dstRef: String?,
    val dstAlt: String?,
    val detail: String,
    val candidateCount: Int,
)

class AppService(
    private val db: Database,
    private val refs: Map<String, RefGenome>,
    private val transcriptsByBuild: Map<String, List<Transcript>>,
    private val transcriptSetVersion: String,
) {
    fun referenceBuilds(): Set<String> = refs.keys

    @Suppress("UNCHECKED_CAST")
    private fun placementFromMap(m: Map<String, Any?>): Placement =
        Placement(m["chrom"] as String, (m["pos"] as Number).toInt(),
            m["ref"] as String, m["alt"] as String, (m["slide"] as Number).toInt())

    private fun placementToMap(p: Placement) = mapOf(
        "chrom" to p.chrom, "pos" to p.pos, "ref" to p.ref, "alt" to p.alt, "slide" to p.slide)

    private fun stepToMap(s: NormStep) = mapOf(
        "kind" to s.kind, "detail" to s.detail, "pos" to s.pos, "ref" to s.ref, "alt" to s.alt)

    @Suppress("UNCHECKED_CAST")
    private fun stepFromMap(m: Map<String, Any?>) = NormStep(
        m["kind"] as String, m["detail"] as String,
        (m["pos"] as Number?)?.toInt(), m["ref"] as String?, m["alt"] as String?)

    private fun sha256(text: String) =
        MessageDigest.getInstance("SHA-256").digest(text.toByteArray())
            .joinToString("") { "%02x".format(it) }

    fun importBatch(batchId: String, build: String, vcf: String): BatchReceipt {
        require(batchId.isNotBlank()) { "batchId is required" }
        val ref = refs[build] ?: error("unknown reference build '$build'; loaded: ${refs.keys}")
        return db.tx {
            findReceipt(batchId)?.let { existing ->
                require(existing.build == build) {
                    "batch '$batchId' already imported against build ${existing.build}, cannot replay under $build"
                }
                return@tx existing.copy(replay = true)
            }
            val lines = VcfParser.parse(vcf)
            var alleles = 0; var mismatches = 0
            val groupIds = LinkedHashSet<Long>()
            lines.forEach { vl ->
                val isComposite = vl.alts.size > 1
                val rawId = insertRaw(vl, batchId, isComposite)
                vl.alts.forEachIndexed { idx, alt ->
                    alleles += 1
                    val v = Variant(vl.chrom, vl.pos, vl.ref, alt)
                    val n = Normalizer.normalize(v, ref)
                    val gid: Long
                    if (n.refMatches) {
                        val c = n.canonical!!
                        gid = upsertGroup(n, ref.build, c, n.placements, n.steps, batchId)
                    } else {
                        mismatches += 1
                        gid = insertMismatchGroup(n, ref.build, batchId)
                    }
                    groupIds += gid
                    addLink(gid, rawId, idx, if (isComposite) "DERIVED" else "ORIGINAL")
                    if (isComposite) addLink(gid, rawId, idx, "ORIGINAL")
                    if (n.refMatches) projectForGroup(gid, ref.build, n.canonical!!, ref)
                }
            }
            val receipt = BatchReceipt(batchId, build, Instant.now().toString(),
                lines.size, alleles, groupIds.size, mismatches, replay = false)
            saveReceipt(receipt, sha256(vcf), vcf)
            receipt
        }
    }

    private fun findReceipt(batchId: String): BatchReceipt? {
        db.conn.prepareStatement(
            "SELECT batch_id,build,received_at,record_count,allele_count,group_count,ref_mismatch_count FROM batches WHERE batch_id=?"
        ).use { ps ->
            ps.setString(1, batchId)
            ps.executeQuery().use { rs ->
                if (!rs.next()) return null
                return BatchReceipt(rs.getString(1), rs.getString(2), rs.getString(3),
                    rs.getInt(4), rs.getInt(5), rs.getInt(6), rs.getInt(7), replay = true)
            }
        }
    }

    private fun saveReceipt(r: BatchReceipt, hash: String, payload: String) {
        db.conn.prepareStatement("""INSERT INTO batches
          (batch_id,build,received_at,record_count,allele_count,group_count,ref_mismatch_count,input_sha256,payload)
          VALUES (?,?,?,?,?,?,?,?,?)""").use { ps ->
            ps.setString(1, r.batchId); ps.setString(2, r.build); ps.setString(3, r.receivedAt)
            ps.setInt(4, r.records); ps.setInt(5, r.alleles); ps.setInt(6, r.groups)
            ps.setInt(7, r.refMismatches); ps.setString(8, hash); ps.setString(9, payload)
            ps.executeUpdate()
        }
    }

    private fun insertRaw(vl: VcfLine, batchId: String, composite: Boolean): Long {
        db.conn.prepareStatement("""INSERT INTO raw_records
          (batch_id,source,line_number,chrom,pos,ref,alt_field,info,raw_line,is_composite)
          VALUES (?,?,?,?,?,?,?,?,?,?)""").use { ps ->
            ps.setString(1, batchId); ps.setString(2, vl.source); ps.setInt(3, vl.lineNumber)
            ps.setString(4, vl.chrom); ps.setInt(5, vl.pos); ps.setString(6, vl.ref)
            ps.setString(7, vl.alts.joinToString(",")); ps.setString(8, vl.info)
            ps.setString(9, vl.rawLine); ps.setInt(10, if (composite) 1 else 0)
            ps.executeUpdate()
            return ps.generatedKeys.use { if (it.next()) it.getLong(1) else error("no raw id") }
        }
    }

    private fun addLink(groupId: Long, rawId: Long, altIndex: Int, role: String) {
        db.conn.prepareStatement("""INSERT OR IGNORE INTO record_links
          (group_id,raw_record_id,alt_index,role) VALUES (?,?,?,?)""").use { ps ->
            ps.setLong(1, groupId); ps.setLong(2, rawId); ps.setInt(3, altIndex)
            ps.setString(4, role); ps.executeUpdate()
        }
    }

    private fun upsertGroup(
        n: NormalizedAllele, build: String, c: Placement,
        placements: List<Placement>, steps: List<NormStep>, batchId: String,
    ): Long {
        val pJson = Json.stringify(placements.map { placementToMap(it) })
        val sJson = Json.stringify(steps.map { stepToMap(it) })
        db.conn.prepareStatement("""INSERT INTO eq_groups
          (build,chrom,canon_pos,canon_ref,canon_alt,expected_ref,
           placement_count,placements_json,norm_steps_json,preference_rule,rules_version,created_batch)
          VALUES (?,?,?,?,?,?,?,?,?,?,?,?)
          ON CONFLICT(build,chrom,canon_pos,canon_ref,canon_alt) DO UPDATE SET
            placement_count=excluded.placement_count,
            placements_json=excluded.placements_json,
            norm_steps_json=excluded.norm_steps_json
          """).use { ps ->
            ps.setString(1, build); ps.setString(2, c.chrom); ps.setInt(3, c.pos)
            ps.setString(4, c.ref); ps.setString(5, c.alt)
            ps.setString(6, n.expectedRef); ps.setInt(7, placements.size)
            ps.setString(8, pJson); ps.setString(9, sJson)
            ps.setString(10, n.preferenceRule); ps.setString(11, ANNOTATION_RULES_VERSION)
            ps.setString(12, batchId) // 12 binds created_batch; literal 0 supplies ref_mismatch
            ps.executeUpdate()
        }
        return groupIdOf(build, c)
    }

    private fun insertMismatchGroup(n: NormalizedAllele, build: String, batchId: String): Long {
        val sJson = Json.stringify(n.steps.map { stepToMap(it) })
        val cols = listOf("build","chrom","canon_pos","canon_ref","canon_alt","ref_mismatch",
            "expected_ref","placement_count","placements_json","norm_steps_json",
            "preference_rule","rules_version","created_batch")
        val marks = List(cols.size) { "?" }.joinToString(",")
        val sql = "INSERT INTO eq_groups (${cols.joinToString(",")}) VALUES ($marks)"
        db.conn.prepareStatement(sql).use { ps ->
            ps.setString(1, build)
            ps.setString(2, n.input.chrom)
            ps.setNull(3, java.sql.Types.INTEGER)
            ps.setNull(4, java.sql.Types.VARCHAR)
            ps.setNull(5, java.sql.Types.VARCHAR)
            ps.setInt(6, 1)
            ps.setString(7, n.expectedRef)
            ps.setInt(8, 0)
            ps.setString(9, "[]")
            ps.setString(10, sJson)
            ps.setString(11, n.preferenceRule)
            ps.setString(12, ANNOTATION_RULES_VERSION)
            ps.setString(13, batchId)
            ps.executeUpdate()
            return ps.generatedKeys.use { if (it.next()) it.getLong(1) else error("no gid") }
        }
    }

    private fun groupIdOf(build: String, c: Placement): Long {
        db.conn.prepareStatement("""SELECT id FROM eq_groups
          WHERE build=? AND chrom=? AND canon_pos=? AND canon_ref=? AND canon_alt=?""").use { ps ->
            ps.setString(1, build); ps.setString(2, c.chrom); ps.setInt(3, c.pos)
            ps.setString(4, c.ref); ps.setString(5, c.alt)
            ps.executeQuery().use { rs ->
                require(rs.next()) { "group upsert lookup failed" }
                return rs.getLong(1)
            }
        }
    }

    private fun projectForGroup(groupId: Long, build: String, canonical: Placement, ref: RefGenome) {
        val txs = transcriptsByBuild[build].orEmpty()
        txs.forEach { tx ->
            val r = TranscriptProjector.project(tx, build, canonical, ref)
            val pathJson = Json.stringify(r.path.map { mapOf("step" to it.step, "detail" to it.detail) })
            db.conn.prepareStatement("""INSERT INTO consequences
              (group_id,tx_id,tx_build,status,failure_type,region,transcript_pos,ref_cdna,alt_cdna,hgvsg_like,path_json)
              VALUES (?,?,?,?,?,?,?,?,?,?,?)
              ON CONFLICT(group_id,tx_id) DO UPDATE SET
                status=excluded.status, failure_type=excluded.failure_type, region=excluded.region,
                transcript_pos=excluded.transcript_pos, ref_cdna=excluded.ref_cdna, alt_cdna=excluded.alt_cdna,
                hgvsg_like=excluded.hgvsg_like, path_json=excluded.path_json""").use { ps ->
                ps.setLong(1, groupId); ps.setString(2, tx.txId); ps.setString(3, tx.build)
                ps.setString(4, r.status.name)
                ps.setString(5, r.failureType); ps.setString(6, r.region?.name)
                if (r.transcriptPos == null) ps.setObject(7, null) else ps.setInt(7, r.transcriptPos)
                ps.setString(8, r.refCdna); ps.setString(9, r.altCdna); ps.setString(10, r.hgvsgLike)
                ps.setString(11, pathJson); ps.executeUpdate()
            }
        }
    }

    data class GroupSummary(
        val id: Long, val build: String, val chrom: String,
        val canonPos: Int?, val canonRef: String?, val canonAlt: String?,
        val refMismatch: Boolean, val placementCount: Int,
        val reviewState: String, val evidenceCount: Int,
    )

    fun listGroups(build: String?): List<GroupSummary> {
        val tail = if (build != null) " WHERE g.build=? ORDER BY g.build,g.chrom,COALESCE(g.canon_pos,0)"
                   else " ORDER BY g.build,g.chrom,COALESCE(g.canon_pos,0)"
        val sql = """SELECT g.id,g.build,g.chrom,g.canon_pos,g.canon_ref,g.canon_alt,
              g.ref_mismatch,g.placement_count,
              COALESCE((SELECT state FROM reviews WHERE group_id=g.id ORDER BY version DESC LIMIT 1),'NONE') state,
              (SELECT COUNT(DISTINCT raw_record_id) FROM record_links WHERE group_id=g.id) ev
            FROM eq_groups g""" + tail
        db.conn.prepareStatement(sql).use { ps ->
            if (build != null) ps.setString(1, build)
            ps.executeQuery().use { rs ->
                val out = mutableListOf<GroupSummary>()
                while (rs.next()) {
                    out += GroupSummary(rs.getLong(1), rs.getString(2), rs.getString(3),
                        rs.getObject(4) as Int?, rs.getString(5), rs.getString(6),
                        rs.getInt(7) == 1, rs.getInt(8), rs.getString(9), rs.getInt(10))
                }
                return out
            }
        }
    }

    private fun consequencesOf(groupId: Long): List<ConsequenceRow> {
        db.conn.prepareStatement("""SELECT tx_id,tx_build,status,failure_type,region,transcript_pos,
              ref_cdna,alt_cdna,hgvsg_like,path_json FROM consequences WHERE group_id=? ORDER BY tx_id""")
            
            .use { ps ->
                ps.setLong(1, groupId)
                ps.executeQuery().use { rs ->
                    val out = mutableListOf<ConsequenceRow>()
                    while (rs.next()) out += ConsequenceRow(
                        rs.getString(1), rs.getString(2), rs.getString(3), rs.getString(4),
                        rs.getString(5), rs.getObject(6) as Int?, rs.getString(7),
                        rs.getString(8), rs.getString(9), rs.getString(10))
                    return out
                }
            }
    }

    private fun evidenceOf(groupId: Long): List<EvidenceRow> {
        db.conn.prepareStatement("""SELECT DISTINCT r.id,r.batch_id,r.source,r.line_number,r.chrom,r.pos,
              r.ref,r.alt_field,l.alt_index,l.role,r.raw_line,r.info
            FROM record_links l JOIN raw_records r ON r.id=l.raw_record_id
            WHERE l.group_id=?
            ORDER BY r.batch_id,r.line_number,l.alt_index,l.role""").use { ps ->
            ps.setLong(1, groupId)
            ps.executeQuery().use { rs ->
                val out = mutableListOf<EvidenceRow>()
                while (rs.next()) out += EvidenceRow(rs.getLong(1), rs.getString(2), rs.getString(3),
                    rs.getInt(4), rs.getString(5), rs.getInt(6), rs.getString(7), rs.getString(8),
                    rs.getInt(9), rs.getString(10), rs.getString(11), rs.getString(12))
                return out
            }
        }
    }

    private fun liftoversOf(groupId: Long): List<LiftoverRowDto> {
        db.conn.prepareStatement("""SELECT to_build,status,dst_chrom,dst_pos,dst_ref,dst_alt,detail,candidate_count
            FROM liftover_results WHERE group_id=? ORDER BY to_build""").use { ps ->
            ps.setLong(1, groupId)
            ps.executeQuery().use { rs ->
                val out = mutableListOf<LiftoverRowDto>()
                while (rs.next()) out += LiftoverRowDto(rs.getString(1), rs.getString(2),
                    rs.getString(3), rs.getObject(4) as Int?, rs.getString(5), rs.getString(6),
                    rs.getString(7), rs.getInt(8))
                return out
            }
        }
    }

    fun getGroup(id: Long): GroupDetail? {
        db.conn.prepareStatement("""SELECT id,build,chrom,canon_pos,canon_ref,canon_alt,ref_mismatch,
              expected_ref,placements_json,norm_steps_json,preference_rule,rules_version
            FROM eq_groups WHERE id=?""").use { ps ->
            ps.setLong(1, id)
            ps.executeQuery().use { rs ->
                if (!rs.next()) return null
                @Suppress("UNCHECKED_CAST")
                val placements = (Json.parse(rs.getString(9)) as List<Any?>).map {
                    placementFromMap(it as Map<String, Any?>)
                }
                @Suppress("UNCHECKED_CAST")
                val steps = (Json.parse(rs.getString(10)) as List<Any?>).map {
                    stepFromMap(it as Map<String, Any?>)
                }
                val (state, version) = latestReview(id)
                return GroupDetail(rs.getLong(1), rs.getString(2), rs.getString(3),
                    rs.getObject(4) as Int?, rs.getString(5), rs.getString(6),
                    rs.getInt(7) == 1, rs.getString(8), placements, steps,
                    rs.getString(11), rs.getString(12),
                    state == "FROZEN", state, evidenceOf(id), consequencesOf(id),
                    liftoversOf(id), version)
            }
        }
    }

    private fun latestReview(groupId: Long): Pair<String, Long> {
        db.conn.prepareStatement("SELECT state,version FROM reviews WHERE group_id=? ORDER BY version DESC LIMIT 1")
            .use { ps ->
                ps.setLong(1, groupId)
                ps.executeQuery().use { rs ->
                    return if (rs.next()) rs.getString(1) to rs.getLong(2) else "NONE" to 0L
                }
            }
    }

    /**
     * Optimistic review submission. [expectedVersion] is the last version the
     * reviewer saw; concurrent commits by another reviewer are rejected.
     */
    fun submitReview(groupId: Long, expectedVersion: Long, reviewer: String, state: String, note: String): Long {
        require(reviewer.isNotBlank()) { "reviewer is required" }
        require(state in listOf("FROZEN", "REF_FLAGGED")) { "state must be FROZEN or REF_FLAGGED" }
        return db.tx {
            val (currentState, currentVersion) = latestReview(groupId)
            require(currentVersion == expectedVersion) {
                "group $groupId changed at version $currentVersion while you reviewed v$expectedVersion; reload and re-submit"
            }
            val next = currentVersion + 1
            db.conn.prepareStatement("""INSERT INTO reviews(group_id,version,reviewer,state,note,submitted_at)
                VALUES (?,?,?,?,?,?)""").use { ps ->
                ps.setLong(1, groupId); ps.setLong(2, next); ps.setString(3, reviewer)
                ps.setString(4, state); ps.setString(5, note); ps.setString(6, Instant.now().toString())
                ps.executeUpdate()
            }
            next
        }
    }

    fun importChain(content: String): ImportImpact = db.tx {
        val (pair, parsed) = ChainIo.parse(content)
        val (fromBuild, toBuild) = pair
        require(fromBuild in refs) { "chain source build '$fromBuild' is not loaded" }
        val warnings = mutableListOf<String>()

        parsed.groupingBy { listOf(it.chrom, it.srcStart, it.srcEnd) }.eachCount()
            .filter { it.value > 1 }.keys.forEach { warnings += "duplicate source interval $it; first segment kept" }

        db.conn.prepareStatement("DELETE FROM liftover_segments WHERE from_build=? AND to_build=?")
            .use { ps -> ps.setString(1, fromBuild); ps.setString(2, toBuild); ps.executeUpdate() }
        db.conn.prepareStatement("DELETE FROM liftover_results WHERE to_build=?").use { ps ->
            ps.setString(1, toBuild); ps.executeUpdate()
        }
        val insertedIds = HashMap<Triple<String, Int, Int>, Long>()
        parsed.forEach { seg ->
            val key = Triple(seg.chrom, seg.srcStart, seg.srcEnd)
            if (key in insertedIds) return@forEach
            db.conn.prepareStatement("""INSERT INTO liftover_segments
              (from_build,to_build,chrom,src_start,src_end,strand,dst_chrom,dst_start,dst_end)
              VALUES (?,?,?,?,?,?,?,?,?)""").use { ps ->
                ps.setString(1, fromBuild); ps.setString(2, toBuild); ps.setString(3, seg.chrom)
                ps.setInt(4, seg.srcStart); ps.setInt(5, seg.srcEnd); ps.setString(6, seg.strand.toString())
                ps.setString(7, seg.dstChrom); ps.setInt(8, seg.dstStart); ps.setInt(9, seg.dstEnd)
                ps.executeUpdate()
                insertedIds[key] = ps.generatedKeys.use { if (it.next()) it.getLong(1) else -1L }
            }
        }

        val stored = parsed.filter { Triple(it.chrom, it.srcStart, it.srcEnd) in insertedIds }
            .map { it.copy(id = insertedIds[Triple(it.chrom, it.srcStart, it.srcEnd)]) }

        var unique = 0; var flipped = 0; var many = 0; var interrupted = 0; var unmapped = 0
        val groups = allGroupPlacements(fromBuild)
        groups.forEach { (groupId, placement) ->
            val r = LiftoverEngine.map(placement, stored)
            when (r.status) {
                LiftoverStatus.UNIQUE -> unique++
                LiftoverStatus.FLIPPED -> flipped++
                LiftoverStatus.MANY_TO_ONE -> many++
                LiftoverStatus.INTERRUPTED -> interrupted++
                LiftoverStatus.UNMAPPED -> unmapped++
            }
            db.conn.prepareStatement("""INSERT INTO liftover_results
              (group_id,to_build,status,dst_chrom,dst_pos,dst_ref,dst_alt,detail,candidate_count)
              VALUES (?,?,?,?,?,?,?,?,?)
              ON CONFLICT(group_id,to_build) DO UPDATE SET
                status=excluded.status,dst_chrom=excluded.dst_chrom,dst_pos=excluded.dst_pos,
                dst_ref=excluded.dst_ref,dst_alt=excluded.dst_alt,detail=excluded.detail,
                candidate_count=excluded.candidate_count""").use { ps ->
                ps.setLong(1, groupId); ps.setString(2, toBuild); ps.setString(3, r.status.name)
                ps.setString(4, r.chrom)
                if (r.pos == null) ps.setObject(5, null) else ps.setInt(5, r.pos)
                ps.setString(6, r.ref); ps.setString(7, r.alt); ps.setString(8, r.detail)
                ps.setInt(9, r.candidateCount); ps.executeUpdate()
            }
        }
        ImportImpact(fromBuild, toBuild, stored.size, unique, flipped, many, interrupted, unmapped, warnings)
    }

    private fun allGroupPlacements(build: String): List<Pair<Long, Placement>> {
        db.conn.prepareStatement("""SELECT id,build,chrom,canon_pos,canon_ref,canon_alt FROM eq_groups
            WHERE build=? AND ref_mismatch=0""").use { ps ->
            ps.setString(1, build)
            ps.executeQuery().use { rs ->
                val out = mutableListOf<Pair<Long, Placement>>()
                while (rs.next()) {
                    out += rs.getLong(1) to Placement(rs.getString(3), rs.getInt(4),
                        rs.getString(5), rs.getString(6), 0)
                }
                return out
            }
        }
    }

    /** Alignment window for the UI: reference context around a placement with a small flank. */
    fun alignmentWindow(build: String, chrom: String, pos: Int, refLen: Int, flank: Int = 12): Map<String, Any?> {
        val genome = refs[build] ?: error("build $build not loaded")
        val start = maxOf(1, pos - flank)
        val end = minOf(genome.length(chrom), pos + refLen - 1 + flank)
        val seq = genome.substring(chrom, start, end - start + 1)
        return mapOf(
            "build" to build, "chrom" to chrom, "start" to start, "end" to end,
            "window" to seq, "alleleStart" to (pos - start), "alleleLength" to refLen,
        )
    }

    fun transcripts(): List<Map<String, Any?>> =
        transcriptsByBuild.flatMap { (build, txs) ->
            txs.map {
                mapOf(
                    "txId" to it.txId, "gene" to it.gene, "chrom" to it.chrom,
                    "strand" to it.strand.toString(), "build" to build,
                    "exons" to it.exons.map { (a, b) -> mapOf("start" to a, "end" to b) },
                    "cdsStart" to it.cdsStart, "cdsEnd" to it.cdsEnd,
                )
            }
        }

    fun batches(): List<Map<String, Any?>> {
        db.conn.prepareStatement("""SELECT batch_id,build,received_at,record_count,allele_count,
              group_count,ref_mismatch_count FROM batches ORDER BY received_at""").use { st ->
            st.executeQuery().use { rs ->
                val out = mutableListOf<Map<String, Any?>>()
                while (rs.next()) {
                    out += mapOf(
                        "batchId" to rs.getString(1), "build" to rs.getString(2),
                        "receivedAt" to rs.getString(3), "records" to rs.getInt(4),
                        "alleles" to rs.getInt(5), "groups" to rs.getInt(6),
                        "refMismatches" to rs.getInt(7),
                    )
                }
                return out
            }
        }
    }
}
