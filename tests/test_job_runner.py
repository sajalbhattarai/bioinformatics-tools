"""
Unit tests for bioinformatics_tools.api.services.job_runner's log-parsing
regexes -- no SSH/SLURM involved, pure string matching.
"""
from bioinformatics_tools.api.services.job_runner import (
    SLURM_SUBMIT_RE, RULE_NAME_FROM_LOG_PATH_RE, SEQUENTIAL_GENOME_RE)


class TestSlurmSubmitRegex:
    """SLURM_SUBMIT_RE only extracts (slurm_id, log_path) now -- genome
    attribution needs the full path to read that job's own log file later
    (ssh_slurm.get_job_genome), so rule-name resolution from the path is a
    separate regex (RULE_NAME_FROM_LOG_PATH_RE) applied to group(2), mirroring
    the actual two-step extraction in run_ssh_task."""

    def test_captures_slurm_id_and_full_log_path(self):
        line = (
            "SLURM jobid 39603778 (log: /path/.snakemake/slurm_logs/"
            "rule_run_quast_batch/39603778.log)."
        )
        slurm_id, log_path = SLURM_SUBMIT_RE.search(line).groups()
        assert slurm_id == "39603778"
        assert log_path == "/path/.snakemake/slurm_logs/rule_run_quast_batch/39603778.log"

    def test_ungrouped_rule_captures_full_rule_name(self):
        log_path = "/path/.snakemake/slurm_logs/rule_run_quast_batch/39603778.log"
        rule_name, group_name = RULE_NAME_FROM_LOG_PATH_RE.search(log_path).groups()
        assert rule_name == "run_quast_batch"
        assert group_name is None

    def test_grouped_rule_captures_only_the_short_group_name(self):
        """run_<tool> and load_<tool>_to_db share one SLURM submission per
        genome (see margie_sb.smk's group: directives) -- the displayed
        name must be the group's own short name ("rasttk"), not the
        snakemake-generated concatenation of every rule in the group
        ("load_rasttk_to_db_run_rasttk"), which read as if only the load
        step ran when the run step's real work happens in the same job."""
        log_path = "/path/.snakemake/slurm_logs/group_rasttk_load_rasttk_to_db_run_rasttk/39608085.log"
        rule_name, group_name = RULE_NAME_FROM_LOG_PATH_RE.search(log_path).groups()
        assert rule_name is None
        assert group_name == "rasttk"

    def test_grouped_rule_kegg(self):
        log_path = "/path/.snakemake/slurm_logs/group_kegg_load_kegg_to_db_run_kegg/39607504.log"
        rule_name, group_name = RULE_NAME_FROM_LOG_PATH_RE.search(log_path).groups()
        assert rule_name is None
        assert group_name == "kegg"

    def test_displayed_name_resolution_prefers_rule_over_group(self):
        """Mirrors the `ungrouped_rule_name or group_name` expression used
        at the actual call site in run_ssh_task."""
        ungrouped_path = "/p/.snakemake/slurm_logs/rule_run_quast_batch/1.log"
        grouped_path = "/p/.snakemake/slurm_logs/group_rasttk_load_rasttk_to_db_run_rasttk/2.log"

        rule_name, group_name = RULE_NAME_FROM_LOG_PATH_RE.search(ungrouped_path).groups()
        assert (rule_name or group_name) == "run_quast_batch"

        rule_name, group_name = RULE_NAME_FROM_LOG_PATH_RE.search(grouped_path).groups()
        assert (rule_name or group_name) == "rasttk"

    def test_real_full_line_with_uuid_jobid_prefix(self):
        """The orchestrator's own internal jobid is a UUID for group jobs
        (GroupJob.jobid), not a small int -- confirmed directly against a
        real dane-api log. SLURM_SUBMIT_RE must not assume it's numeric."""
        line = (
            "[2026-06-23 04:53:32] INFO bioinformatics_tools.workflow_tools.workflow: "
            "[snakemake] Job f2bd3300-d2be-5a07-91a9-b299cb448c4f has been submitted "
            "with SLURM jobid 39684417 (log: /scratch/negishi/bhattar3/margie/output/"
            "2026-06-23-0453/.snakemake/slurm_logs/"
            "group_rasttk_load_rasttk_to_db_run_rasttk/39684417.log)."
        )
        slurm_id, log_path = SLURM_SUBMIT_RE.search(line).groups()
        assert slurm_id == "39684417"
        rule_name, group_name = RULE_NAME_FROM_LOG_PATH_RE.search(log_path).groups()
        assert (rule_name or group_name) == "rasttk"


class TestSequentialGenomeRegex:
    """margie_sb's sequential per-organism orchestrator (workflow.py's
    _run_pipeline_batch_sequential) prints one marker line per genome
    transition -- see this module's own SEQUENTIAL_GENOME_RE comment for
    why that's needed alongside STEPS_PROGRESS_RE."""

    def test_matches_marker_line(self):
        line = "=== SEQUENTIAL: genome 3/12 (Genus_species) phase4-8 starting ==="
        match = SEQUENTIAL_GENOME_RE.search(line)
        assert match is not None
        assert match.groups() == ("3", "12")

    def test_no_match_on_unrelated_line(self):
        line = "5 of 20 steps (25%) done"
        assert SEQUENTIAL_GENOME_RE.search(line) is None
