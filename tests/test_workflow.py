"""
Tests for WorkflowBase: _run_subprocess, build_executable,
_parse_snakemake_output, _run_pipeline, do_quick_example, do_fresh_test.

All tests are mocked — no snakemake installation required.
"""
import subprocess
import threading
from contextlib import contextmanager
from unittest.mock import MagicMock, patch

import pytest

from bioinformatics_tools.workflow_tools import workflow as workflow_module
from bioinformatics_tools.workflow_tools.workflow import WorkflowBase, WORKFLOWS, _scripts_versioned
from bioinformatics_tools.workflow_tools.workflow_registry import WORKFLOWS


# ---------------------------------------------------------------------------
# Fixture: create a WorkflowBase without triggering CLI __init__
# ---------------------------------------------------------------------------

@pytest.fixture
def wf():
    """Build a WorkflowBase instance without CLI init."""
    obj = WorkflowBase.__new__(WorkflowBase)
    conf = MagicMock()
    conf.get = MagicMock(side_effect=lambda key, default=None: {
        'main_database': '/tmp/test-margie.db',
    }.get(key, default))
    obj.conf = conf
    obj.report = None
    obj.workflow_id = 'test'
    obj.timestamp = '010101-0000'
    return obj


# ---------------------------------------------------------------------------
# _parse_snakemake_output
# ---------------------------------------------------------------------------

class TestParseSnakemakeOutput:

    def test_parses_failed_rules(self):
        stderr = (
            "Error in rule run_pfam:\n"
            "    some details\n"
            "Error in rule run_cog:\n"
            "    more details\n"
            "2 of 5 steps (40%) done\n"
        )
        result = WorkflowBase._parse_snakemake_output(stderr)
        assert result['completed'] == 2
        assert result['total'] == 5
        assert result['failed'] == 2
        assert set(result['failed_rules']) == {'run_pfam', 'run_cog'}

    def test_full_success(self):
        stderr = "5 of 5 steps (100%) done\n"
        result = WorkflowBase._parse_snakemake_output(stderr)
        assert result['completed'] == 5
        assert result['total'] == 5
        assert result['failed'] == 0
        assert result['failed_rules'] == []

    def test_empty_stderr(self):
        result = WorkflowBase._parse_snakemake_output("")
        assert result == {'total': 0, 'completed': 0, 'failed': 0, 'failed_rules': []}

    def test_failed_rules_without_steps_line(self):
        stderr = "Error in rule step_flaky:\n    jobid: 2\n"
        result = WorkflowBase._parse_snakemake_output(stderr)
        assert result['failed_rules'] == ['step_flaky']
        assert result['failed'] == 1
        # total estimated from completed (0) + failed (1)
        assert result['total'] == 1


# ---------------------------------------------------------------------------
# _run_subprocess
# ---------------------------------------------------------------------------

class TestRunSubprocess:

    def test_success_returns_completed_process(self, wf):
        mock_proc = MagicMock()
        mock_proc.returncode = 0
        mock_proc.stdout = iter(['ok\n'])
        mock_proc.stderr = iter([])
        mock_proc.wait.return_value = None

        with patch('bioinformatics_tools.workflow_tools.workflow.subprocess.Popen', return_value=mock_proc):
            result = wf._run_subprocess(['snakemake', '-s', 'test.smk'])
        assert result.returncode == 0
        assert result.stdout == 'ok'

    def test_nonzero_still_returns(self, wf):
        mock_proc = MagicMock()
        mock_proc.returncode = 1
        mock_proc.stdout = iter([])
        mock_proc.stderr = iter(['Error in rule x:\n'])
        mock_proc.wait.return_value = None

        with patch('bioinformatics_tools.workflow_tools.workflow.subprocess.Popen', return_value=mock_proc):
            result = wf._run_subprocess(['snakemake', '-s', 'test.smk'])
        assert result.returncode == 1

    def test_launch_failure_returns_none(self, wf):
        with patch('bioinformatics_tools.workflow_tools.workflow.subprocess.Popen',
                   side_effect=FileNotFoundError('snakemake not found')):
            result = wf._run_subprocess(['snakemake', '-s', 'test.smk'])
        assert result is None
        # failed() should have been called
        assert wf.report is not None
        assert wf.report.status.indicates_failure


# ---------------------------------------------------------------------------
# build_executable
# ---------------------------------------------------------------------------

class TestBuildExecutable:

    def test_has_keep_going(self, wf):
        key = WORKFLOWS['selftest']
        cmd = wf.build_executable(key, mode='dev')
        assert '--keep-going' in cmd

    def test_dev_mode_no_slurm_executor(self, wf):
        key = WORKFLOWS['selftest']
        cmd = wf.build_executable(key, mode='dev')
        assert '--executor=slurm' not in cmd

    def test_non_dev_has_slurm_executor(self, wf):
        key = WORKFLOWS['selftest']
        cmd = wf.build_executable(key, mode='notdev')
        # Should appear exactly once
        assert cmd.count('--executor=slurm') == 1

    def test_config_dict_appended(self, wf):
        key = WORKFLOWS['selftest']
        cmd = wf.build_executable(key, config_overrides={'foo': 'bar', 'baz': '42'}, mode='dev')
        assert '--config' in cmd
        idx = cmd.index('--config')
        assert 'foo=bar' in cmd[idx + 1:]
        assert 'baz=42' in cmd[idx + 1:]

    def test_dev_mode_no_default_resources(self, wf):
        key = WORKFLOWS['selftest']
        cmd = wf.build_executable(key, mode='dev')
        assert '--default-resources' not in cmd

    def test_no_target_leaves_command_unchanged(self, wf):
        """Checks that no target adds neither --nolock nor a target token."""
        key = WORKFLOWS['selftest']
        cmd = wf.build_executable(key, mode='dev')
        assert '--nolock' not in cmd
        assert 'rasttk_all' not in cmd
        assert 'phase4_8_one_genome' not in cmd

    def test_target_appears_before_config(self, wf):
        """Checks that the target precedes --config, which consumes every
        following token as key=value."""
        key = WORKFLOWS['selftest']
        cmd = wf.build_executable(key, config_overrides={'foo': 'bar'}, mode='dev', target='rasttk_all')
        assert 'rasttk_all' in cmd
        assert cmd.index('rasttk_all') < cmd.index('--config')

    def test_target_adds_nolock(self, wf):
        key = WORKFLOWS['selftest']
        cmd = wf.build_executable(key, mode='dev', target='phase4_8_one_genome')
        assert '--nolock' in cmd
        assert cmd.index('--nolock') < cmd.index('phase4_8_one_genome')


# ---- _start_background_subprocess ----

class TestStartBackgroundSubprocess:

    def test_returns_immediately_without_waiting(self, wf):
        mock_proc = MagicMock()
        mock_proc.stdout = iter([])
        mock_proc.stderr = iter([])
        mock_proc.poll.return_value = None  # still "running"

        with patch('bioinformatics_tools.workflow_tools.workflow.subprocess.Popen', return_value=mock_proc):
            handle = wf._start_background_subprocess(['snakemake', '-s', 'test.smk'])

        assert handle is not None
        assert handle.poll() is None
        mock_proc.wait.assert_not_called()

    def test_exposes_real_returncode_once_finished(self, wf):
        mock_proc = MagicMock()
        mock_proc.stdout = iter([])
        mock_proc.stderr = iter([])
        mock_proc.poll.return_value = 0
        mock_proc.wait.return_value = 0

        with patch('bioinformatics_tools.workflow_tools.workflow.subprocess.Popen', return_value=mock_proc):
            handle = wf._start_background_subprocess(['snakemake', '-s', 'test.smk'])

        assert handle.wait() == 0
        assert handle.poll() == 0

    def test_launch_failure_returns_none(self, wf):
        with patch('bioinformatics_tools.workflow_tools.workflow.subprocess.Popen',
                   side_effect=FileNotFoundError('snakemake not found')):
            handle = wf._start_background_subprocess(['snakemake', '-s', 'test.smk'])
        assert handle is None


# ---------------------------------------------------------------------------
# _run_pipeline
# ---------------------------------------------------------------------------

class TestRunPipeline:

    def test_unknown_key_fails(self, wf):
        ret = wf._run_pipeline('nonexistent_workflow', {})
        assert ret == 1
        assert wf.report is not None
        assert wf.report.status.indicates_failure

    @patch('bioinformatics_tools.workflow_tools.workflow.cache_sif_files')
    def test_success_path(self, mock_cache, wf):
        fake_proc = subprocess.CompletedProcess(
            args=['snakemake'], returncode=0,
            stdout='Building DAG\n', stderr='5 of 5 steps (100%) done\n',
        )
        with patch.object(wf, '_run_subprocess', return_value=fake_proc):
            wf._run_pipeline('example', {'input_fasta': 'test.fa'})

        assert wf.report is not None
        assert wf.report.status.indicates_success
        assert wf.report.data['workflow'] == 'example'
        assert wf.report.data['returncode'] == 0
        assert wf.report.data['rules_summary']['completed'] == 5

    @patch('bioinformatics_tools.workflow_tools.workflow.cache_sif_files')
    def test_failure_path_does_not_call_succeeded(self, mock_cache, wf):
        """Regression test: when snakemake fails, self.succeeded() must NOT be called."""
        fake_proc = subprocess.CompletedProcess(
            args=['snakemake'], returncode=1,
            stdout='', stderr='Error in rule run_pfam:\n    jobid: 3\n1 of 3 steps (33%) done\n',
        )
        with patch.object(wf, '_run_subprocess', return_value=fake_proc):
            ret = wf._run_pipeline('example', {'input_fasta': 'test.fa'})

        assert ret == 1
        assert wf.report.status.indicates_failure
        assert wf.report.data['rules_summary']['failed_rules'] == ['run_pfam']

    @patch('bioinformatics_tools.workflow_tools.workflow.cache_sif_files')
    def test_launch_failure_returns_early(self, mock_cache, wf):
        with patch.object(wf, '_run_subprocess', return_value=None):
            ret = wf._run_pipeline('example', {'input_fasta': 'test.fa'})
        assert ret == 1

    @patch('bioinformatics_tools.workflow_tools.workflow.cache_sif_files')
    @patch('bioinformatics_tools.workflow_tools.workflow.log_workflow_run')
    @patch('bioinformatics_tools.workflow_tools.workflow.store_all')
    def test_store_all_skipped_on_failure(self, mock_store, mock_log, mock_cache, wf):
        fake_proc = subprocess.CompletedProcess(
            args=['snakemake'], returncode=1, stdout='', stderr='Error in rule x:\n',
        )
        cache_map = {'prodigal': ['out.tkn']}
        smk_config = {'input_fasta': 'test.fa', 'main_database': '/tmp/test.db'}
        with patch.object(wf, '_run_subprocess', return_value=fake_proc):
            wf._run_pipeline('example', smk_config, cache_map)
        mock_store.assert_not_called()
        mock_log.assert_called_once()
        assert mock_log.call_args.kwargs['status'] == 'failed'

    @patch('bioinformatics_tools.workflow_tools.workflow.cache_sif_files')
    @patch('bioinformatics_tools.workflow_tools.workflow.log_workflow_run')
    @patch('bioinformatics_tools.workflow_tools.workflow.store_all')
    @patch('bioinformatics_tools.workflow_tools.workflow.restore_all', return_value={})
    def test_store_all_called_on_success(self, mock_restore, mock_store, mock_log, mock_cache, wf):
        fake_proc = subprocess.CompletedProcess(
            args=['snakemake'], returncode=0, stdout='', stderr='',
        )
        cache_map = {'prodigal': ['out.tkn']}
        smk_config = {'input_fasta': 'test.fa', 'main_database': '/tmp/test.db'}
        with patch.object(wf, '_run_subprocess', return_value=fake_proc):
            wf._run_pipeline('example', smk_config, cache_map)
        mock_store.assert_called_once()
        mock_log.assert_called_once()
        assert mock_log.call_args.kwargs['status'] == 'success'

    def test_selftest_skips_cache_sif(self, wf):
        """selftest has empty sif_files, so cache_sif_files should not be called."""
        fake_proc = subprocess.CompletedProcess(
            args=['snakemake'], returncode=0, stdout='', stderr='',
        )
        with patch('bioinformatics_tools.workflow_tools.workflow.cache_sif_files') as mock_cache, \
             patch.object(wf, '_run_subprocess', return_value=fake_proc):
            wf._run_pipeline('selftest', {'workdir': '/tmp'}, mode='dev')
        mock_cache.assert_not_called()

    @patch('bioinformatics_tools.workflow_tools.workflow.cache_sif_files',
           side_effect=__import__('bioinformatics_tools.workflow_tools.bapptainer',
                                  fromlist=['CacheSifError']).CacheSifError('download failed'))
    def test_cache_sif_failure(self, mock_cache, wf):
        ret = wf._run_pipeline('example', {'input_fasta': 'test.fa'})
        assert ret == 1
        assert wf.report.status.indicates_failure

    @patch('bioinformatics_tools.workflow_tools.workflow.log_workflow_run')
    @patch('bioinformatics_tools.workflow_tools.workflow.store_all')
    @patch('bioinformatics_tools.workflow_tools.workflow.restore_all', return_value={})
    def test_pipeline_with_input_file_key(self, mock_restore, mock_store, mock_log, wf):
        """_run_pipeline uses input_file key when input_fasta is absent (selftest path)."""
        fake_proc = subprocess.CompletedProcess(
            args=['snakemake'], returncode=0, stdout='', stderr='',
        )
        cache_map = {'step_a': ['step_a/sample-a-step_a.out']}
        smk_config = {'input_file': '/tmp/sample-a.txt', 'main_database': '/tmp/sample.db'}
        with patch.object(wf, '_run_subprocess', return_value=fake_proc):
            wf._run_pipeline('selftest', smk_config, cache_map, mode='dev')
        mock_restore.assert_called_once_with('/tmp/sample.db', '/tmp/sample-a.txt', cache_map)
        mock_log.assert_called_once()
        assert mock_log.call_args.kwargs['status'] == 'success'


# ---- _run_pipeline_batch_sequential ----
# Path and time.sleep are mocked so the polling loops run instantly, and every
# process the runner starts is replaced by a finished-looking fake (see
# _fake_launchers), so no snakemake is needed.

def _target_genome_of(cmd: list[str]) -> str | None:
    """Returns target_genome=<x> from a built snakemake command's --config section."""
    for tok in cmd:
        if tok.startswith('target_genome='):
            return tok.split('=', 1)[1]
    return None


class _FakeBackground:
    """Stands in for workflow._BackgroundProcess (.poll(), .wait(),
    .proc.returncode, .stdout_lines, .stderr_lines).

    Finished as soon as it is created, or, given `running_until`, running until
    that Event is set (so Stage 1 can be kept alive while tokens appear)."""

    def __init__(self, returncode: int = 0, stderr: str = '',
                 running_until: threading.Event | None = None):
        self.proc = MagicMock(returncode=returncode)
        self.stdout_lines: list[str] = []
        self.stderr_lines = stderr.splitlines()
        self._returncode = returncode
        self._running_until = running_until

    def poll(self):
        if self._running_until is not None and not self._running_until.is_set():
            return None
        return self._returncode

    def wait(self):
        if self._running_until is not None:
            self._running_until.wait()
        return self._returncode


@contextmanager
def _fake_launchers(wf, *, stage2_returncodes=(), stage1_running_until=None):
    """Replaces every process _run_pipeline_batch_sequential starts and records
    the commands it was asked to run:

      rec['stage1']  Stage 1 wave launches (_start_background_subprocess, no target_genome)
      rec['stage2']  per-genome Stage 2 launches (target_genome=<name>), in start order
      rec['run']     _run_subprocess calls (Stage 3 llm_all and the finalize targets)

    Stage 2 launch n exits with stage2_returncodes[n] (0 when not given); a
    non-zero one carries a snakemake-style "Error in rule" stderr. Stage 1
    waves exit at once unless stage1_running_until is given. The real
    _Stage1Waves drives the waves; it is only subclassed to record instances so
    its thread can be joined before the patches are undone."""
    rec = {'stage1': [], 'stage2': [], 'run': []}

    def fake_start_background(cmd):
        if _target_genome_of(cmd) is None:
            rec['stage1'].append(cmd)
            return _FakeBackground(running_until=stage1_running_until)
        index = len(rec['stage2'])
        rec['stage2'].append(cmd)
        rc = stage2_returncodes[index] if index < len(stage2_returncodes) else 0
        return _FakeBackground(returncode=rc, stderr='Error in rule run_pfam:\n' if rc else '')

    def fake_run_subprocess(cmd):
        rec['run'].append(cmd)
        return subprocess.CompletedProcess(args=cmd, returncode=0, stdout='', stderr='')

    waves = []

    class _RecordingWaves(workflow_module._Stage1Waves):
        def __init__(self, *args, **kwargs):
            waves.append(self)
            super().__init__(*args, **kwargs)

    with patch.object(wf, '_start_background_subprocess', side_effect=fake_start_background), \
         patch.object(wf, '_run_subprocess', side_effect=fake_run_subprocess), \
         patch('bioinformatics_tools.workflow_tools.workflow._Stage1Waves', _RecordingWaves), \
         patch('bioinformatics_tools.workflow_tools.workflow.subprocess.run') as mock_run:
        # The finalize reorganize step calls subprocess.run directly.
        mock_run.return_value = subprocess.CompletedProcess(args=[], returncode=0, stdout='', stderr='')
        try:
            yield rec
        finally:
            if stage1_running_until is not None:
                stage1_running_until.set()
            for wave in waves:
                wave.wait()


class TestRunPipelineBatchSequential:
    """Stage 1 runs as waves of background snakemake processes over explicit
    file targets (_Stage1Waves, --jobs=1); Stage 2 is one background process
    per genome (rule phase4_12_one_genome with target_genome=<name>) started in
    RASTtk-token order; Stage 3 (llm_all) and the finalize targets run through
    _run_subprocess afterwards."""

    SMK_CONFIG = {'output_dir': '/tmp/seqtest', 'main_database': ''}

    @patch('bioinformatics_tools.workflow_tools.workflow.cache_sif_files')
    @patch('bioinformatics_tools.workflow_tools.workflow.locate_local_sif_files')
    @patch('bioinformatics_tools.workflow_tools.workflow.time.sleep')
    def test_stage1_targets_rasttk_all_stage2_targets_phase4_8_one_genome(self, mock_sleep, mock_locate, mock_cache, wf):
        """Checks that Stage 1 and Stage 2 request the right margie_sb.smk targets:
        Stage 1 the genome's RASTtk/GTDB-Tk db tokens as file targets with
        --jobs=1, Stage 2 rule phase4_12_one_genome with target_genome=<name>,
        then the two finalize targets (no llm_all while LLM is off)."""
        genome_files = {'genome1': '/g1.fa'}

        with patch('bioinformatics_tools.workflow_tools.workflow.Path') as mock_path_cls, \
             _fake_launchers(wf) as rec:
            mock_path_cls.return_value.exists.return_value = True  # GTDB-Tk outputs already on disk
            wf._run_pipeline_batch_sequential('margie_sb', dict(self.SMK_CONFIG), genome_files,
                                              cache_map_fn=lambda s: {}, mode='dev')

        # Stage 1: GTDB-Tk is already available, so a single RASTtk wave over file targets.
        assert len(rec['stage1']) == 1
        stage1_cmd = rec['stage1'][0]
        assert '--jobs=1' in stage1_cmd
        assert '--nolock' in stage1_cmd
        assert any(tok.endswith('/genome1/rasttk/rasttk_db.tkn') for tok in stage1_cmd)
        assert any(tok.endswith('/genome1/gtdbtk/gtdbtk_db.tkn') for tok in stage1_cmd)
        assert 'rasttk_all' not in stage1_cmd

        # Stage 2: one process for the one genome.
        assert len(rec['stage2']) == 1
        assert 'phase4_12_one_genome' in rec['stage2'][0]
        assert _target_genome_of(rec['stage2'][0]) == 'genome1'

        # LLM off: no Stage 3; the finalize targets follow Stage 2.
        assert not any('llm_all' in cmd for cmd in rec['run'])
        assert sum('queue_sqlite_backup_snapshot' in cmd for cmd in rec['run']) == 1
        assert sum('run_report_figures_global' in cmd for cmd in rec['run']) == 1
        assert wf.report.status.indicates_success

    @patch('bioinformatics_tools.workflow_tools.workflow.cache_sif_files')
    @patch('bioinformatics_tools.workflow_tools.workflow.locate_local_sif_files')
    @patch('bioinformatics_tools.workflow_tools.workflow.time.sleep')
    def test_llm_on_stops_stage2_before_llm_and_runs_llm_all(self, mock_sleep, mock_locate, mock_cache, wf):
        """Checks that with run_llm on, Stage 2 targets phase4_12_one_genome_no_llm
        and one Stage 3 llm_all run follows."""
        genome_files = {'genome1': '/g1.fa'}

        with patch('bioinformatics_tools.workflow_tools.workflow.Path') as mock_path_cls, \
             _fake_launchers(wf) as rec:
            mock_path_cls.return_value.exists.return_value = True
            wf._run_pipeline_batch_sequential('margie_sb', {**self.SMK_CONFIG, 'run_llm': True}, genome_files,
                                              cache_map_fn=lambda s: {}, mode='dev')

        assert len(rec['stage2']) == 1
        assert 'phase4_12_one_genome_no_llm' in rec['stage2'][0]
        assert sum('llm_all' in cmd for cmd in rec['run']) == 1
        assert wf.report.status.indicates_success

    @patch('bioinformatics_tools.workflow_tools.workflow.cache_sif_files')
    @patch('bioinformatics_tools.workflow_tools.workflow.locate_local_sif_files')
    @patch('bioinformatics_tools.workflow_tools.workflow.time.sleep')
    def test_genomes_processed_in_token_detection_order(self, mock_sleep, mock_locate, mock_cache, wf):
        """Checks that genomes enter Stage 2 in the order their RASTtk token
        appears; time.sleep's side_effect separates poll rounds, and Stage 1
        stays alive until the last token is there."""
        genome_files = {'genome1': '/g1.fa', 'genome2': '/g2.fa', 'genome3': '/g3.fa'}

        ready = {'genome2'}
        stage1_done = threading.Event()

        def fake_sleep(seconds):
            if 'genome1' not in ready:
                ready.add('genome1')
            elif 'genome3' not in ready:
                ready.add('genome3')
                stage1_done.set()  # Stage 1 exits once the last token exists
        mock_sleep.side_effect = fake_sleep

        def fake_path(path_str):
            m = MagicMock()
            m.exists.return_value = any(g in str(path_str) for g in ready)
            return m

        with patch('bioinformatics_tools.workflow_tools.workflow.Path', side_effect=fake_path), \
             _fake_launchers(wf, stage1_running_until=stage1_done) as rec:
            wf._run_pipeline_batch_sequential('margie_sb', dict(self.SMK_CONFIG), genome_files,
                                              cache_map_fn=lambda s: {}, mode='dev')

        assert [_target_genome_of(cmd) for cmd in rec['stage2']] == ['genome2', 'genome1', 'genome3']
        # Three waves: the GTDB-Tk batch, RASTtk for genome2 (already had
        # GTDB-Tk outputs), then RASTtk for genome1/genome3 after the batch.
        assert len(rec['stage1']) == 3
        assert wf.report.status.indicates_success

    @patch('bioinformatics_tools.workflow_tools.workflow.cache_sif_files')
    @patch('bioinformatics_tools.workflow_tools.workflow.locate_local_sif_files')
    @patch('bioinformatics_tools.workflow_tools.workflow.time.sleep')
    def test_failure_halts_before_next_genome(self, mock_sleep, mock_locate, mock_cache, wf):
        """Checks that a Stage 2 failure halts the queue: no third genome is
        started, nothing runs after it, and the run is reported as failed."""
        genome_files = {'genome1': '/g1.fa', 'genome2': '/g2.fa', 'genome3': '/g3.fa'}

        with patch('bioinformatics_tools.workflow_tools.workflow.Path') as mock_path_cls, \
             _fake_launchers(wf, stage2_returncodes=(0, 1)) as rec:
            mock_path_cls.return_value.exists.return_value = True  # every genome is RASTtk-ready
            ret = wf._run_pipeline_batch_sequential('margie_sb', dict(self.SMK_CONFIG), genome_files,
                                                     cache_map_fn=lambda s: {}, mode='dev')

        assert len(rec['stage2']) == 2
        assert rec['run'] == []  # no Stage 3 / finalize after a halt
        assert ret == 1
        assert wf.report.status.indicates_failure

    @patch('bioinformatics_tools.workflow_tools.workflow.cache_sif_files')
    @patch('bioinformatics_tools.workflow_tools.workflow.locate_local_sif_files')
    @patch('bioinformatics_tools.workflow_tools.workflow.time.sleep')
    def test_genome_that_never_becomes_ready_is_skipped_not_blocking(self, mock_sleep, mock_locate, mock_cache, wf):
        """Checks that a genome whose RASTtk token never appears is skipped once
        Stage 1 has exited, without blocking the others or failing the run."""
        genome_files = {'genome1': '/g1.fa', 'genome2': '/g2.fa'}

        def fake_path(path_str):
            m = MagicMock()
            m.exists.return_value = 'genome2' in str(path_str)  # genome1 never becomes ready
            return m

        with patch('bioinformatics_tools.workflow_tools.workflow.Path', side_effect=fake_path), \
             _fake_launchers(wf) as rec:
            ret = wf._run_pipeline_batch_sequential('margie_sb', dict(self.SMK_CONFIG), genome_files,
                                                     cache_map_fn=lambda s: {}, mode='dev')

        assert [_target_genome_of(cmd) for cmd in rec['stage2']] == ['genome2']
        assert ret is None  # finished(), not failed()
        assert wf.report.status.indicates_success
        assert not wf.report.status.indicates_failure
        assert wf.report.data['genomes_processed'] == 1
        assert wf.report.data['genomes_skipped'] == ['genome1']


# ---- do_margie_sb ----

class TestDoMargieSb:

    def test_genome_cache_map_covers_every_phase4_8_tool(self, wf):
        """Checks that do_margie_sb passes a per-genome cache_map_fn covering every
        phase4-8 tool, with the full file set for multi-file tools, plus the
        script-versioned consolidation/labeling steps."""
        wf.conf.get = MagicMock(side_effect=lambda key, default=None: {
            'input': '/genomes/genome1.fa',
            'main_database': '/tmp/test-margie.db',
        }.get(key, default))

        captured = {}

        def fake_run_pipeline_batch_sequential(key_name, smk_config, genome_files, **kwargs):
            captured['cache_map_fn'] = kwargs['cache_map_fn']

        entitlement = {'usage_type': 'academic', 'licensed_tools': [], 'source': 'test'}
        with patch('bioinformatics_tools.workflow_tools.workflow.discover_genomes',
                   return_value={'genome1': '/genomes/genome1.fa'}), \
             patch('bioinformatics_tools.workflow_tools.env_keys.ensure_api_keys'), \
             patch('bioinformatics_tools.workflow_tools.license_gate.ensure_cli_license',
                   return_value=entitlement), \
             patch.object(wf, '_run_pipeline_batch_sequential', side_effect=fake_run_pipeline_batch_sequential):
            wf.do_margie_sb(mode='dev')

        cache_map = captured['cache_map_fn']('genome1')

        consolidation = _scripts_versioned('consolidation')
        labeling = _scripts_versioned('labeling')
        assert consolidation.startswith('consolidation@')
        assert labeling.startswith('labeling@')
        expected_tools = {
            'cog', 'pfam', 'merops', 'tcdb', 'uniprot', 'kegg', 'eggnog', 'dbcan',
            'pgap', 'geneprop', 'operon', 'tmbed', 'signalp6', 'deepsig', 'psortb',
            'signalp4', 'tigrfam', 'phobius', 'envelope', 'interpro',
            consolidation, labeling,
        }
        assert set(cache_map.keys()) == expected_tools
        # quast/gtdbtk/rasttk are excluded (see _genome_cache_map).
        assert 'quast' not in cache_map
        assert 'gtdbtk' not in cache_map
        assert 'rasttk' not in cache_map

        assert len(cache_map['cog']) == 2  # results.tsv + token
        assert len(cache_map['tigrfam']) == 3  # results.tsv + domtbl.out + token
        assert len(cache_map['phobius']) == 3  # results.tsv + top1.tsv + token
        assert len(cache_map['envelope']) == 3  # results.tsv + summary.tsv + token
        assert len(cache_map['interpro']) == 38  # (results.tsv + token) + 18 dbs * 2 files each
        assert len(cache_map[consolidation]) == 6  # tables + manifest + compute/db tokens
        assert len(cache_map[labeling]) == 6

        for paths in cache_map.values():
            for path in paths:
                assert path.endswith(('.tsv', '.tkn', '.out', '.json'))


# ---------------------------------------------------------------------------
# do_quick_example
# ---------------------------------------------------------------------------

class TestDoQuickExample:

    @patch('bioinformatics_tools.workflow_tools.workflow.log_workflow_run')
    @patch('bioinformatics_tools.workflow_tools.workflow.store_all')
    @patch('bioinformatics_tools.workflow_tools.workflow.restore_all', return_value={})
    def test_quick_example_passes_cache_map(self, mock_restore, mock_store, mock_log, wf):
        """do_quick_example should call _run_pipeline with a cache_map matching the step keys."""
        fake_proc = subprocess.CompletedProcess(
            args=['snakemake'], returncode=0, stdout='', stderr='',
        )
        with patch.object(wf, '_run_subprocess', return_value=fake_proc):
            wf.do_quick_example()

        # restore_all was called with a cache_map containing all step keys
        call_args = mock_restore.call_args
        cache_map = call_args[0][2]
        assert set(cache_map.keys()) == {
            'step_a', 'step_a_db', 'step_b', 'step_b_db', 'step_c', 'step_c_db',
        }
        # Each value should be a list of output path strings
        assert len(cache_map['step_a']) == 2  # .out and .extra
        assert len(cache_map['step_c']) == 2  # .tsv and _count.tsv

        # store_all and log_workflow_run should both be called on success
        mock_store.assert_called_once()
        mock_log.assert_called_once()
        assert mock_log.call_args.kwargs['status'] == 'success'

    def test_quick_example_uses_selftest_workflow_key(self, wf):
        """do_quick_example runs the 'selftest' workflow key (no sif files)."""
        fake_proc = subprocess.CompletedProcess(
            args=['snakemake'], returncode=0, stdout='', stderr='',
        )
        with patch('bioinformatics_tools.workflow_tools.workflow.restore_all', return_value={}), \
             patch('bioinformatics_tools.workflow_tools.workflow.store_all'), \
             patch.object(wf, '_run_subprocess', return_value=fake_proc) as mock_sub, \
             patch('bioinformatics_tools.workflow_tools.workflow.cache_sif_files') as mock_cache:
            wf.do_quick_example()

        # selftest has no sif_files, so cache_sif_files should not be called
        mock_cache.assert_not_called()
        # _run_subprocess was called (snakemake command was built)
        mock_sub.assert_called_once()


# ---------------------------------------------------------------------------
# do_fresh_test
# ---------------------------------------------------------------------------

class TestDoFreshTest:

    @patch('bioinformatics_tools.workflow_tools.workflow.log_workflow_run')
    @patch('bioinformatics_tools.workflow_tools.workflow.store_all')
    @patch('bioinformatics_tools.workflow_tools.workflow.restore_all', return_value={})
    def test_fresh_test_uses_cache_map(self, mock_restore, mock_store, mock_log, wf):
        """do_fresh_test passes cache_map and uses real margie_db for store/restore."""
        fake_proc = subprocess.CompletedProcess(
            args=['snakemake'], returncode=0, stdout='', stderr='',
        )
        with patch.object(wf, '_run_subprocess', return_value=fake_proc):
            wf.do_fresh_test()

        # restore_all, store_all, and log_workflow_run should all be called on success
        mock_restore.assert_called_once()
        mock_store.assert_called_once()
        mock_log.assert_called_once()
        assert mock_log.call_args.kwargs['status'] == 'success'

        # cache_map should have all step keys
        cache_map = mock_restore.call_args[0][2]
        assert set(cache_map.keys()) == {
            'step_a', 'step_a_db', 'step_b', 'step_b_db', 'step_c', 'step_c_db',
        }

    @patch('bioinformatics_tools.workflow_tools.workflow.store_all')
    @patch('bioinformatics_tools.workflow_tools.workflow.restore_all', return_value={})
    def test_fresh_test_passes_inject_failure(self, mock_restore, mock_store, wf):
        """do_fresh_test should pass inject_failure through to smk_config."""
        fake_proc = subprocess.CompletedProcess(
            args=['snakemake'], returncode=0, stdout='', stderr='',
        )
        with patch.object(wf, '_run_subprocess', return_value=fake_proc) as mock_sub:
            wf.do_fresh_test(inject_failure=True)

        # Check the snakemake command includes inject_failure=true in config
        cmd = mock_sub.call_args[0][0]
        config_str = ' '.join(cmd)
        assert 'inject_failure=true' in config_str

    @patch('bioinformatics_tools.workflow_tools.workflow.store_all')
    @patch('bioinformatics_tools.workflow_tools.workflow.restore_all', return_value={})
    def test_fresh_test_runs_selftest_key(self, mock_restore, mock_store, wf):
        """do_fresh_test uses the 'selftest' workflow key."""
        fake_proc = subprocess.CompletedProcess(
            args=['snakemake'], returncode=0, stdout='', stderr='',
        )
        with patch.object(wf, '_run_subprocess', return_value=fake_proc) as mock_sub:
            wf.do_fresh_test()

        cmd = mock_sub.call_args[0][0]
        # Should reference selftest.smk
        assert any('selftest.smk' in arg for arg in cmd)


# ---------------------------------------------------------------------------
# _build_result
# ---------------------------------------------------------------------------

class TestBuildResult:

    def test_builds_structured_dict(self, wf):
        proc = subprocess.CompletedProcess(
            args=['snakemake'], returncode=1,
            stdout='some output', stderr='Error in rule bad_rule:\n2 of 3 steps (66%) done\n',
        )
        result = wf._build_result('margie', proc)
        assert result['workflow'] == 'margie'
        assert result['returncode'] == 1
        assert result['rules_summary']['failed_rules'] == ['bad_rule']
        assert result['rules_summary']['completed'] == 2
        assert 'some output' in result['stdout_tail']
