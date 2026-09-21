package app

import java.nio.file.Files
import java.nio.file.Path
import java.sql.Connection
import java.sql.Statement

object Schema {
    const val VERSION = 1

    val DDL = """
CREATE TABLE IF NOT EXISTS schema_meta(
  key TEXT PRIMARY KEY,
  value TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS batches(
  batch_id TEXT PRIMARY KEY,
  build TEXT NOT NULL,
  received_at TEXT NOT NULL,
  record_count INTEGER NOT NULL,
  allele_count INTEGER NOT NULL,
  group_count INTEGER NOT NULL,
  ref_mismatch_count INTEGER NOT NULL,
  input_sha256 TEXT NOT NULL,
  payload TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS refs(
  build TEXT PRIMARY KEY,
  fasta_path TEXT NOT NULL,
  loaded_at TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS raw_records(
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  batch_id TEXT NOT NULL REFERENCES batches(batch_id),
  source TEXT NOT NULL,
  line_number INTEGER NOT NULL,
  chrom TEXT NOT NULL,
  pos INTEGER NOT NULL,
  ref TEXT NOT NULL,
  alt_field TEXT NOT NULL,
  info TEXT NOT NULL,
  raw_line TEXT NOT NULL,
  is_composite INTEGER NOT NULL,
  UNIQUE(batch_id, line_number)
);

CREATE TABLE IF NOT EXISTS eq_groups(
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  build TEXT NOT NULL,
  chrom TEXT NOT NULL,
  canon_pos INTEGER,
  canon_ref TEXT,
  canon_alt TEXT,
  ref_mismatch INTEGER NOT NULL DEFAULT 0,
  expected_ref TEXT,
  placement_count INTEGER NOT NULL DEFAULT 1,
  placements_json TEXT NOT NULL,
  norm_steps_json TEXT NOT NULL,
  preference_rule TEXT NOT NULL,
  rules_version TEXT NOT NULL,
  created_batch TEXT NOT NULL,
  UNIQUE(build, chrom, canon_pos, canon_ref, canon_alt)
);

CREATE TABLE IF NOT EXISTS record_links(
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  group_id INTEGER NOT NULL REFERENCES eq_groups(id),
  raw_record_id INTEGER NOT NULL REFERENCES raw_records(id),
  alt_index INTEGER NOT NULL,
  role TEXT NOT NULL CHECK(role IN ('ORIGINAL','DERIVED'))
);

CREATE TABLE IF NOT EXISTS consequences(
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  group_id INTEGER NOT NULL REFERENCES eq_groups(id),
  tx_id TEXT NOT NULL,
  tx_build TEXT NOT NULL,
  status TEXT NOT NULL,
  failure_type TEXT,
  region TEXT,
  transcript_pos INTEGER,
  ref_cdna TEXT,
  alt_cdna TEXT,
  hgvsg_like TEXT,
  path_json TEXT NOT NULL,
  UNIQUE(group_id, tx_id)
);

CREATE TABLE IF NOT EXISTS liftover_segments(
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  from_build TEXT NOT NULL,
  to_build TEXT NOT NULL,
  chrom TEXT NOT NULL,
  src_start INTEGER NOT NULL,
  src_end INTEGER NOT NULL,
  strand TEXT NOT NULL,
  dst_chrom TEXT NOT NULL,
  dst_start INTEGER NOT NULL,
  dst_end INTEGER NOT NULL,
  UNIQUE(from_build, to_build, chrom, src_start, src_end)
);

CREATE TABLE IF NOT EXISTS liftover_results(
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  group_id INTEGER NOT NULL REFERENCES eq_groups(id),
  to_build TEXT NOT NULL,
  status TEXT NOT NULL,
  dst_chrom TEXT,
  dst_pos INTEGER,
  dst_ref TEXT,
  dst_alt TEXT,
  detail TEXT NOT NULL,
  candidate_count INTEGER NOT NULL,
  UNIQUE(group_id, to_build)
);

CREATE TABLE IF NOT EXISTS reviews(
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  group_id INTEGER NOT NULL REFERENCES eq_groups(id),
  version INTEGER NOT NULL,
  reviewer TEXT NOT NULL,
  state TEXT NOT NULL CHECK(state IN ('NONE','FROZEN','REF_FLAGGED')),
  note TEXT NOT NULL DEFAULT '',
  submitted_at TEXT NOT NULL,
  UNIQUE(group_id, version)
);

CREATE INDEX IF NOT EXISTS idx_raw_batch ON raw_records(batch_id);
CREATE UNIQUE INDEX IF NOT EXISTS uq_original_evidence
  ON record_links(group_id, raw_record_id) WHERE role='ORIGINAL';
CREATE INDEX IF NOT EXISTS idx_links_raw ON record_links(raw_record_id);
CREATE INDEX IF NOT EXISTS idx_links_group ON record_links(group_id);
CREATE INDEX IF NOT EXISTS idx_cons_group ON consequences(group_id);
CREATE INDEX IF NOT EXISTS idx_groups_build ON eq_groups(build);
"""
}

class Database(val conn: Connection) {
    companion object {
        fun open(path: Path): Database {
            Files.createDirectories(path.toAbsolutePath().parent)
            val c = java.sql.DriverManager.getConnection("jdbc:sqlite:" + path.toAbsolutePath())
            c.autoCommit = false
            c.createStatement().use { it.execute("PRAGMA foreign_keys=ON") }
            val db = Database(c)
            db.init()
            return db
        }
    }

    fun init() {
        conn.createStatement().use { st ->
            Schema.DDL.trimIndent().split(";").map { it.trim() }
                .filter { it.isNotEmpty() }
                .forEach { st.execute(it) }
        }
        val existing = conn.prepareStatement("SELECT value FROM schema_meta WHERE key='schema_version'")
            .use { ps -> ps.executeQuery().use { rs -> if (rs.next()) rs.getString(1) else null } }
        if (existing == null) {
            conn.prepareStatement("INSERT INTO schema_meta(key,value) VALUES('schema_version',?)")
                .use { it.setString(1, Schema.VERSION.toString()); it.executeUpdate() }
            conn.commit()
        } else {
            check(existing == Schema.VERSION.toString()) { "unexpected schema version $existing" }
        }
    }

    fun <T> tx(body: () -> T): T {
        val oldAuto = conn.autoCommit
        conn.autoCommit = false
        try {
            val r = body()
            conn.commit()
            return r
        } catch (t: Throwable) {
            conn.rollback()
            throw t
        } finally {
            conn.autoCommit = oldAuto
        }
    }
}
