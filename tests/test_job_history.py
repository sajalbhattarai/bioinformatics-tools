"""
Unit tests for bioinformatics_tools.api.services.job_history, run against a
real temporary SQLite database rather than a mocked sqlite3.
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
        assert "owner_username" in cols
        assert "owner_cluster_username" in cols
        assert "selected_tools" in cols
        assert "relaunched_from" in cols

    def test_pre_existing_old_schema_table_gets_migrated(self, tmp_path):
        """Checks that an older table gains the new columns via ALTER TABLE and keeps its rows."""
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
        assert "owner_username" in cols
        assert "owner_cluster_username" in cols
        assert "selected_tools" in cols
        assert "relaunched_from" in cols

        row = job_history.get_job(db_path, "old-job")
        assert row is not None
        assert row["status"] == "completed"
        assert row["selected_tools"] is None

    def test_ensure_table_idempotent(self, tmp_path):
        """Checks that calling ensure_table twice does not fail."""
        db_path = str(tmp_path / "twice.db")
        job_history.ensure_table(db_path)
        job_history.ensure_table(db_path)
        assert "selected_tools" in _columns(db_path)


class TestRecordJobCreatedRoundTrip:
    def test_selected_tools_and_relaunched_from_persist(self, tmp_path):
        db_path = str(tmp_path / "round_trip.db")
        job_history.record_job_created(
            db_path, "job-1", "margie_sb", "/genomes/e.fasta", "/out/2026-06-21-1200",
            owner_username="alice", owner_cluster_username="alice_cluster",
            selected_tools="quast,gtdbtk", relaunched_from="job-0",
        )
        row = job_history.get_job(db_path, "job-1")
        assert row["owner_username"] == "alice"
        assert row["owner_cluster_username"] == "alice_cluster"
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
    """Tests the logs/slurm_jobs/containers snapshot written by job_store.finalize(),
    which must round-trip as Python lists rather than JSON-encoded strings."""

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

    def test_list_jobs_leaves_out_snapshot_fields_get_job_keeps_them(self, tmp_path):
        """Checks that list_jobs omits the snapshot lists and trims the log, while get_job keeps them."""
        db_path = str(tmp_path / "snapshot_list.db")
        job_history.record_job_created(
            db_path, "job-4", "margie_sb", "/genomes/e.fasta", "/out/2026-06-21-1200",
        )
        snapshot = [{"job_id": "222", "rule": "rasttk", "status": "COMPLETED", "time": "01:00:00"}]
        job_history.record_job_updated(db_path, "job-4", slurm_jobs=snapshot, logs="x" * 20000 + "the end")

        rows = job_history.list_jobs(db_path)
        assert rows[0]["slurm_jobs"] == []
        assert rows[0]["logs"].endswith("the end") and len(rows[0]["logs"]) == job_history.LIST_LOG_TAIL
        row = job_history.get_job(db_path, "job-4")
        assert row["slurm_jobs"] == snapshot and len(row["logs"]) == 20007

    def test_no_snapshot_yet_decodes_to_empty_list(self, tmp_path):
        """Checks that NULL slurm_jobs/containers columns decode to []."""
        db_path = str(tmp_path / "no_snapshot.db")
        job_history.record_job_created(
            db_path, "job-5", "margie_sb", "/genomes/e.fasta", "/out/2026-06-21-1200",
        )
        row = job_history.get_job(db_path, "job-5")
        assert row["slurm_jobs"] == []
        assert row["containers"] == []
        assert row["logs"] is None
        assert row["relaunched_from"] is None


class TestOwnershipFiltering:
    def test_list_and_get_are_scoped_by_owner(self, tmp_path):
        db_path = str(tmp_path / "owners.db")
        job_history.record_job_created(
            db_path, "job-a", "margie_sb", "/genomes/a.fasta", "/scratch/team/alice/out",
            owner_username="alice", owner_cluster_username="alice",
        )
        job_history.record_job_created(
            db_path, "job-b", "margie_sb", "/genomes/b.fasta", "/scratch/team/bob/out",
            owner_username="bob", owner_cluster_username="bob",
        )

        rows = job_history.list_jobs(
            db_path,
            owner_username="alice",
            owner_cluster_username="alice",
        )
        assert [r["job_id"] for r in rows] == ["job-a"]

        assert job_history.get_job(
            db_path,
            "job-b",
            owner_username="alice",
            owner_cluster_username="alice",
        ) is None

    def test_legacy_ownerless_rows_can_fallback_to_cluster_path(self, tmp_path):
        db_path = str(tmp_path / "legacy_ownerless.db")
        job_history.record_job_created(
            db_path, "legacy-a", "margie_sb", "/scratch/cluster/alice/genomes/a.fasta",
            "/scratch/cluster/alice/margie-output/2026-08-04-1200",
        )
        job_history.record_job_created(
            db_path, "legacy-b", "margie_sb", "/scratch/cluster/bob/genomes/b.fasta",
            "/scratch/cluster/bob/margie-output/2026-08-04-1200",
        )

        rows = job_history.list_jobs(
            db_path,
            owner_username="alice-app",
            owner_cluster_username="alice",
        )
        assert [r["job_id"] for r in rows] == ["legacy-a"]
