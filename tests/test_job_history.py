"""
Unit tests for bioinformatics_tools.api.services.job_history.

This module runs directly on the cluster account against a real SQLite
file (see its own module docstring) -- so these tests use a real (temp)
SQLite database rather than mocking sqlite3, unlike the SSH-exec-based
utility tests elsewhere in this suite.
"""
import sqlite3

from bioinformatics_tools.api.services import job_history


def _columns(db_path: str) -> set[str]:
    conn = sqlite3.connect(db_path)
    try:
        return {row[1] for row in conn.execute("PRAGMA table_info(api_jobs)")}
    finally:
        conn.close()


class TestEnsureTableMigration:
    def test_fresh_db_gets_new_columns(self, tmp_path):
        db_path = str(tmp_path / "fresh.db")
        job_history.ensure_table(db_path)
        cols = _columns(db_path)
        assert "selected_tools" in cols
        assert "relaunched_from" in cols

    def test_pre_existing_old_schema_table_gets_migrated(self, tmp_path):
        """A table created before selected_tools/relaunched_from existed
        must gain both columns via ALTER TABLE, without losing existing
        rows -- CREATE TABLE IF NOT EXISTS alone wouldn't touch an
        already-existing table's schema at all."""
        db_path = str(tmp_path / "old.db")
        conn = sqlite3.connect(db_path)
        conn.execute("""
            CREATE TABLE api_jobs (
                job_id TEXT PRIMARY KEY,
                workflow TEXT NOT NULL,
                genome_path TEXT,
                output_dir TEXT,
                work_dir TEXT,
                status TEXT NOT NULL,
                phase TEXT,
                created_at TEXT NOT NULL,
                updated_at TEXT NOT NULL
            )
        """)
        conn.execute(
            "INSERT INTO api_jobs (job_id, workflow, status, created_at, updated_at) "
            "VALUES ('old-job', 'margie_sb', 'completed', 'then', 'then')"
        )
        conn.commit()
        conn.close()

        job_history.ensure_table(db_path)

        cols = _columns(db_path)
        assert "selected_tools" in cols
        assert "relaunched_from" in cols

        row = job_history.get_job(db_path, "old-job")
        assert row is not None
        assert row["status"] == "completed"
        assert row["selected_tools"] is None

    def test_ensure_table_idempotent(self, tmp_path):
        """Calling ensure_table twice (e.g. on every record_job_created
        call, as it already is) must not error on the second pass."""
        db_path = str(tmp_path / "twice.db")
        job_history.ensure_table(db_path)
        job_history.ensure_table(db_path)
        assert "selected_tools" in _columns(db_path)


class TestRecordJobCreatedRoundTrip:
    def test_selected_tools_and_relaunched_from_persist(self, tmp_path):
        db_path = str(tmp_path / "round_trip.db")
        job_history.record_job_created(
            db_path, "job-1", "margie_sb", "/genomes/e.fasta", "/out/2026-06-21-1200",
            selected_tools="quast,gtdbtk", relaunched_from="job-0",
        )
        row = job_history.get_job(db_path, "job-1")
        assert row["selected_tools"] == "quast,gtdbtk"
        assert row["relaunched_from"] == "job-0"

    def test_defaults_to_none_when_not_given(self, tmp_path):
        db_path = str(tmp_path / "defaults.db")
        job_history.record_job_created(
            db_path, "job-2", "margie_sb", "/genomes/e.fasta", "/out/2026-06-21-1200",
        )
        row = job_history.get_job(db_path, "job-2")
        assert row["selected_tools"] is None


class TestFinalSnapshotPersistence:
    """logs/slurm_jobs/containers: job_store.finalize()'s one-time snapshot,
    routed through record_job_updated -- slurm_jobs/containers must
    round-trip as real Python lists, not JSON-encoded strings, since
    job_history.py decodes them before returning (so the SSH/JSON
    transport layer sees proper nested arrays, not double-encoded text)."""

    def test_logs_slurm_jobs_containers_round_trip_as_real_lists(self, tmp_path):
        db_path = str(tmp_path / "snapshot.db")
        job_history.record_job_created(
            db_path, "job-3", "margie_sb", "/genomes/e.fasta", "/out/2026-06-21-1200",
        )
        job_history.record_job_updated(
            db_path, "job-3",
            status="completed", phase="Done",
            logs="line one\nline two\n",
            slurm_jobs=[{"job_id": "111", "rule": "quast", "status": "COMPLETED", "time": "00:01:00"}],
            containers=[{"name": "quast", "version": "5.0"}],
        )

        row = job_history.get_job(db_path, "job-3")
        assert row["logs"] == "line one\nline two\n"
        assert row["slurm_jobs"] == [{"job_id": "111", "rule": "quast", "status": "COMPLETED", "time": "00:01:00"}]
        assert row["containers"] == [{"name": "quast", "version": "5.0"}]

    def test_list_jobs_also_decodes_snapshot_fields(self, tmp_path):
        db_path = str(tmp_path / "snapshot_list.db")
        job_history.record_job_created(
            db_path, "job-4", "margie_sb", "/genomes/e.fasta", "/out/2026-06-21-1200",
        )
        job_history.record_job_updated(
            db_path, "job-4", slurm_jobs=[{"job_id": "222", "rule": "rasttk", "status": "COMPLETED", "time": "01:00:00"}],
        )

        rows = job_history.list_jobs(db_path)
        assert rows[0]["slurm_jobs"] == [{"job_id": "222", "rule": "rasttk", "status": "COMPLETED", "time": "01:00:00"}]

    def test_no_snapshot_yet_decodes_to_empty_list(self, tmp_path):
        """A job that hasn't reached finalize() yet (still running, or
        created before this feature existed) has NULL slurm_jobs/
        containers columns -- must decode to [], not crash on json.loads(None)."""
        db_path = str(tmp_path / "no_snapshot.db")
        job_history.record_job_created(
            db_path, "job-5", "margie_sb", "/genomes/e.fasta", "/out/2026-06-21-1200",
        )
        row = job_history.get_job(db_path, "job-5")
        assert row["slurm_jobs"] == []
        assert row["containers"] == []
        assert row["logs"] is None
        assert row["relaunched_from"] is None
