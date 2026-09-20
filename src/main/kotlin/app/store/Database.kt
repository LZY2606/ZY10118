package app.store

import app.model.*
import java.io.File
import java.sql.Connection
import java.sql.DriverManager
import java.sql.ResultSet

/** SQLite-backed state store.
 *
 * Concurrency/atomicity contract:
 *  - every batch import and every reference/map update runs inside a single
 *    serialized transaction (BEGIN IMMEDIATE); either all rows commit or the
 *    database is left untouched;
 *  - retrying the same batch ID returns the original receipt;
 *  - reviews use an optimistic version column, so a stale review cannot
 *    overwrite a newer submission.
 */
class Database(path: String) {
    val connection: Connection

    init {
        Class.forName("org.sqlite.JDBC")
        File(path).absoluteFile.parentFile?.mkdirs()
        connection = DriverManager.getConnection("jdbc:sqlite:$path")
        connection.createStatement().use { st ->
            st.execute("PRAGMA journal_mode=WAL")
            st.execute("PRAGMA foreign_keys=ON")
            st.execute("PRAGMA busy_timeout=10000")
        }
        migrate()
        connection.autoCommit = false
    }

    private fun migrate() {
        connection.createStatement().use { st ->
            st.executeUpdate(
                """
                CREATE TABLE IF NOT EXISTS builds (
                  name TEXT PRIMARY KEY,
                  imported_at TEXT NOT NULL,
                  fasta_path TEXT NOT NULL
                )""".trimIndent()
            )
            st.executeUpdate(
                """
                CREATE TABLE IF NOT EXISTS contigs (
                  build TEXT NOT NULL,
                  contig TEXT NOT NULL,
                  length INTEGER NOT NULL,
                  sequence TEXT NOT NULL,
                  PRIMARY KEY (build, contig)
                )""".trimIndent()
            )
            st.executeUpdate(
                """
                CREATE TABLE IF NOT EXISTS transcripts (
                  build TEXT NOT NULL,
                  transcript_id TEXT NOT NULL,
                  contig TEXT NOT NULL,
                  strand TEXT NOT NULL,
                  exon_starts TEXT NOT NULL,
                  exon_ends TEXT NOT NULL,
                  cds_start INTEGER NOT NULL,
                  cds_end INTEGER NOT NULL,
                  cds_exon_starts TEXT NOT NULL,
                  cds_exon_ends TEXT NOT NULL,
                  PRIMARY KEY (build, transcript_id)
                )""".trimIndent()
            )
            st.executeUpdate(
                """
                CREATE TABLE IF NOT EXISTS batches (
                  batch_id TEXT PRIMARY KEY,
                  build TEXT NOT NULL,
                  payload_sha256 TEXT NOT NULL,
                  received_at TEXT NOT NULL,
                  record_count INTEGER NOT NULL,
                  receipt TEXT NOT NULL
                )""".trimIndent()
            )
            st.executeUpdate(
                """
                CREATE TABLE IF NOT EXISTS observations (
                  id INTEGER PRIMARY KEY AUTOINCREMENT,
                  batch_id TEXT NOT NULL,
                  line_number INTEGER NOT NULL,
                  contig TEXT NOT NULL,
                  pos INTEGER NOT NULL,
                  ref TEXT NOT NULL,
                  alt TEXT NOT NULL,
                  build TEXT NOT NULL,
                  ref_status TEXT NOT NULL,
                  parent_line_number INTEGER,
                  alt_index INTEGER NOT NULL,
                  child_count INTEGER NOT NULL,
                  source_id TEXT NOT NULL,
                  canonical_id TEXT,
                  UNIQUE (batch_id, line_number, alt_index)
                )""".trimIndent()
            )
            st.executeUpdate(
                """
                CREATE TABLE IF NOT EXISTS canonical_groups (
                  canonical_id TEXT PRIMARY KEY,
                  build TEXT NOT NULL,
                  contig TEXT NOT NULL,
                  pos INTEGER NOT NULL,
                  ref TEXT NOT NULL,
                  alt TEXT NOT NULL,
                  rules_version TEXT NOT NULL,
                  equivalent_anchors TEXT NOT NULL,
                  unique_placement INTEGER NOT NULL,
                  norm_steps TEXT NOT NULL,
                  frozen INTEGER NOT NULL DEFAULT 0,
                  frozen_anchor INTEGER,
                  marked_mismatch INTEGER NOT NULL DEFAULT 0
                )""".trimIndent()
            )
            st.executeUpdate(
                """
                CREATE TABLE IF NOT EXISTS projections (
                  id INTEGER PRIMARY KEY AUTOINCREMENT,
                  canonical_id TEXT NOT NULL,
                  transcript_id TEXT NOT NULL,
                  build TEXT NOT NULL,
                  rules_version TEXT NOT NULL,
                  consequence TEXT NOT NULL,
                  failure TEXT,
                  cds_start INTEGER,
                  cds_end INTEGER,
                  protein_start INTEGER,
                  protein_end INTEGER,
                  ref_codon TEXT,
                  alt_codon TEXT,
                  ref_aa TEXT,
                  alt_aa TEXT,
                  path_json TEXT NOT NULL,
                  affected_exons TEXT NOT NULL,
                  UNIQUE (canonical_id, transcript_id)
                )""".trimIndent()
            )
            st.executeUpdate(
                """
                CREATE TABLE IF NOT EXISTS liftovers (
                  id INTEGER PRIMARY KEY AUTOINCREMENT,
                  source_build TEXT NOT NULL,
                  dest_build TEXT NOT NULL,
                  imported_at TEXT NOT NULL,
                  payload_sha256 TEXT NOT NULL,
                  UNIQUE (source_build, dest_build)
                )""".trimIndent()
            )
            st.executeUpdate(
                """
                CREATE TABLE IF NOT EXISTS liftover_segments (
                  id INTEGER PRIMARY KEY AUTOINCREMENT,
                  liftover_id INTEGER NOT NULL,
                  ord INTEGER NOT NULL,
                  src_contig TEXT NOT NULL,
                  src_start INTEGER NOT NULL,
                  src_end INTEGER NOT NULL,
                  dst_contig TEXT NOT NULL,
                  dst_start INTEGER NOT NULL,
                  dst_end INTEGER NOT NULL,
                  dst_strand TEXT NOT NULL
                )""".trimIndent()
            )
            st.executeUpdate(
                """
                CREATE TABLE IF NOT EXISTS impact_reports (
                  id INTEGER PRIMARY KEY AUTOINCREMENT,
                  liftover_id INTEGER NOT NULL,
                  canonical_id TEXT NOT NULL,
                  status TEXT NOT NULL,
                  dest_contig TEXT,
                  dest_pos INTEGER,
                  dest_ref TEXT,
                  dest_alt TEXT,
                  segments_used TEXT NOT NULL,
                  note TEXT NOT NULL,
                  created_at TEXT NOT NULL
                )""".trimIndent()
            )
            st.executeUpdate(
                """
                CREATE TABLE IF NOT EXISTS reviews (
                  id INTEGER PRIMARY KEY AUTOINCREMENT,
                  canonical_id TEXT NOT NULL,
                  version INTEGER NOT NULL,
                  reviewer TEXT NOT NULL,
                  action TEXT NOT NULL,
                  anchor_override INTEGER,
                  note TEXT,
                  submitted_at TEXT NOT NULL,
                  UNIQUE (canonical_id, version)
                )""".trimIndent()
            )
        }
    }

    fun <T> transaction(block: () -> T): T {
        // Keep one connection with autoCommit disabled. The first write opens
        // a SQLite transaction; busy_timeout plus WAL serialize writers, and
        // commit/rollback make each import/update atomic.
        connection.autoCommit = false
        try {
            val result = block()
            connection.commit()
            return result
        } catch (t: Throwable) {
            connection.rollback()
            throw t
        }
    }

    fun close() = connection.close()
}

fun ResultSet.getStringOrNull(col: String): String? = getString(col)
