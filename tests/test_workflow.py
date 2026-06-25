"""
Tests for WorkflowBase: _run_subprocess, build_executable,
_parse_snakemake_output, _run_pipeline, do_quick_example, do_fresh_test.

All tests are mocked — no snakemake installation required.
"""
import subprocess
from unittest.mock import MagicMock, patch

import pytest

from bioinformatics_tools.workflow_tools.workflow import WorkflowBase, WORKFLOWS
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
        """No target (the default) must not add --nolock or any target token --
        every existing caller's command stays byte-for-byte the same."""
        key = WORKFLOWS['selftest']
        cmd = wf.build_executable(key, mode='dev')
        assert '--nolock' not in cmd
        assert 'rasttk_all' not in cmd
        assert 'phase4_8_one_genome' not in cmd

    def test_target_appears_before_config(self, wf):
        """--config greedily consumes every following token as a key=value
        pair (confirmed directly against a real snakemake invocation) --
        the target must come before it or snakemake crashes on a bare
        rule name."""
        key = WORKFLOWS['selftest']
        cmd = wf.build_executable(key, config_overrides={'foo': 'bar'}, mode='dev', target='rasttk_all')
        assert 'rasttk_all' in cmd
        assert cmd.index('rasttk_all') < cmd.index('--config')

    def test_target_adds_nolock(self, wf):
        key = WORKFLOWS['selftest']
        cmd = wf.build_executable(key, mode='dev', target='phase4_8_one_genome')
        assert '--nolock' in cmd
        assert cmd.index('--nolock') < cmd.index('phase4_8_one_genome')


# ---------------------------------------------------------------------------
# _start_background_subprocess
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# _run_pipeline_batch_sequential
#
# margie_sb's sequential per-organism orchestrator: Stage 1 (rasttk_all,
# all genomes, backgrounded) + Stage 2 (phase4_8_one_genome, one genome at
# a time, in the order each genome's RASTtk token actually appears).
# Path and time.sleep are patched at module level (not real-filesystem /
# real-wait) so these run instantly and deterministically; _run_subprocess
# and _start_background_subprocess are the same seams TestRunSubprocess/
# TestStartBackgroundSubprocess exercise directly, mocked here instead so
# the orchestration loop itself is what's under test.
# ---------------------------------------------------------------------------

def _target_genome_of(cmd: list[str]) -> str | None:
    """Pull target_genome=<x> out of a built snakemake command's --config
    section, the same shape build_executable() produces."""
    for tok in cmd:
        if tok.startswith('target_genome='):
            return tok.split('=', 1)[1]
    return None


class TestRunPipelineBatchSequential:

    @staticmethod
    def _success_proc():
        return subprocess.CompletedProcess(args=['snakemake'], returncode=0, stdout='', stderr='')

    @staticmethod
    def _failure_proc():
        return subprocess.CompletedProcess(
            args=['snakemake'], returncode=1, stdout='', stderr='Error in rule run_pfam:\n')

    @patch('bioinformatics_tools.workflow_tools.workflow.cache_sif_files')
    @patch('bioinformatics_tools.workflow_tools.workflow.locate_local_sif_files')
    @patch('bioinformatics_tools.workflow_tools.workflow.time.sleep')
    def test_stage1_targets_rasttk_all_stage2_targets_phase4_8_one_genome(self, mock_sleep, mock_locate, mock_cache, wf):
        """Smoke test for the actual rule names threaded through -- the
        ordering/failure tests below mock these calls away entirely, so
        this is the one place that checks Stage 1/Stage 2 ask for the
        right margie_sb.smk targets."""
        genome_files = {'genome1': '/g1.fa'}
        smk_config = {'output_dir': '/tmp/seqtest', 'main_database': ''}

        fake_background = MagicMock()
        fake_background.poll.return_value = None
        background_calls = []

        def fake_start_background(cmd):
            background_calls.append(cmd)
            return fake_background

        stage2_calls = []

        def fake_run_subprocess(cmd):
            stage2_calls.append(cmd)
            return self._success_proc()

        with patch('bioinformatics_tools.workflow_tools.workflow.Path') as mock_path_cls, \
             patch.object(wf, '_start_background_subprocess', side_effect=fake_start_background), \
             patch.object(wf, '_run_subprocess', side_effect=fake_run_subprocess):
            mock_path_cls.return_value.exists.return_value = True
            wf._run_pipeline_batch_sequential('margie_sb', smk_config, genome_files,
                                              cache_map_fn=lambda s: {}, mode='dev')

        assert len(background_calls) == 1
        assert 'rasttk_all' in background_calls[0]
        assert len(stage2_calls) == 1
        assert 'phase4_8_one_genome' in stage2_calls[0]
        assert _target_genome_of(stage2_calls[0]) == 'genome1'
        assert wf.report.status.indicates_success

    @patch('bioinformatics_tools.workflow_tools.workflow.cache_sif_files')
    @patch('bioinformatics_tools.workflow_tools.workflow.locate_local_sif_files')
    @patch('bioinformatics_tools.workflow_tools.workflow.time.sleep')
    def test_genomes_processed_in_token_detection_order(self, mock_sleep, mock_locate, mock_cache, wf):
        """Genomes must enter Stage 2 in the order their RASTtk token
        actually appears -- not genome_files' dict order, not alphabetical.
        genome2 is made "ready" immediately; genome1 and genome3 only
        become ready on later poll rounds (simulated via time.sleep's
        side_effect, since it's the only thing separating poll rounds)."""
        genome_files = {'genome1': '/g1.fa', 'genome2': '/g2.fa', 'genome3': '/g3.fa'}
        smk_config = {'output_dir': '/tmp/seqtest', 'main_database': ''}

        ready = {'genome2'}

        def fake_sleep(seconds):
            if 'genome1' not in ready:
                ready.add('genome1')
            elif 'genome3' not in ready:
                ready.add('genome3')
        mock_sleep.side_effect = fake_sleep

        def fake_path(path_str):
            m = MagicMock()
            m.exists.return_value = any(g in str(path_str) for g in ready)
            return m

        fake_background = MagicMock()
        fake_background.poll.return_value = None

        processed_order = []

        def fake_run_subprocess(cmd):
            processed_order.append(_target_genome_of(cmd))
            return self._success_proc()

        with patch('bioinformatics_tools.workflow_tools.workflow.Path', side_effect=fake_path), \
             patch.object(wf, '_start_background_subprocess', return_value=fake_background), \
             patch.object(wf, '_run_subprocess', side_effect=fake_run_subprocess):
            wf._run_pipeline_batch_sequential('margie_sb', smk_config, genome_files,
                                              cache_map_fn=lambda s: {}, mode='dev')

        assert processed_order == ['genome2', 'genome1', 'genome3']
        assert wf.report.status.indicates_success

    @patch('bioinformatics_tools.workflow_tools.workflow.cache_sif_files')
    @patch('bioinformatics_tools.workflow_tools.workflow.locate_local_sif_files')
    @patch('bioinformatics_tools.workflow_tools.workflow.time.sleep')
    def test_failure_halts_before_next_genome(self, mock_sleep, mock_locate, mock_cache, wf):
        """A real Stage 2 failure must halt the queue -- later genomes,
        even ones already RASTtk-ready, must never be attempted. All three
        genomes are made "ready" immediately (set iteration order across
        pending decides which two of the three get processed before the
        failure on the 2nd call -- doesn't matter which two; what matters
        is there's never a 3rd call)."""
        genome_files = {'genome1': '/g1.fa', 'genome2': '/g2.fa', 'genome3': '/g3.fa'}
        smk_config = {'output_dir': '/tmp/seqtest', 'main_database': ''}

        fake_background = MagicMock()
        fake_background.poll.return_value = None

        call_count = {'n': 0}

        def fake_run_subprocess(cmd):
            call_count['n'] += 1
            return self._failure_proc() if call_count['n'] == 2 else self._success_proc()

        with patch('bioinformatics_tools.workflow_tools.workflow.Path') as mock_path_cls, \
             patch.object(wf, '_start_background_subprocess', return_value=fake_background), \
             patch.object(wf, '_run_subprocess', side_effect=fake_run_subprocess):
            mock_path_cls.return_value.exists.return_value = True
            ret = wf._run_pipeline_batch_sequential('margie_sb', smk_config, genome_files,
                                                     cache_map_fn=lambda s: {}, mode='dev')

        assert call_count['n'] == 2
        assert ret == 1
        assert wf.report.status.indicates_failure

    @patch('bioinformatics_tools.workflow_tools.workflow.cache_sif_files')
    @patch('bioinformatics_tools.workflow_tools.workflow.locate_local_sif_files')
    @patch('bioinformatics_tools.workflow_tools.workflow.time.sleep')
    def test_genome_that_never_becomes_ready_is_skipped_not_blocking(self, mock_sleep, mock_locate, mock_cache, wf):
        """If Stage 1 exits without ever producing genome1's RASTtk token
        (a real phase1-3 failure for that genome specifically -- --keep-going
        already let Stage 1 continue past it), genome1 must be skipped, not
        block genome2 from being processed, and not fail the whole run."""
        genome_files = {'genome1': '/g1.fa', 'genome2': '/g2.fa'}
        smk_config = {'output_dir': '/tmp/seqtest', 'main_database': ''}

        poll_calls = {'n': 0}

        def fake_poll():
            poll_calls['n'] += 1
            return None if poll_calls['n'] == 1 else 0  # still running once, then exited

        fake_background = MagicMock()
        fake_background.poll.side_effect = fake_poll
        fake_background.wait.return_value = 0

        def fake_path(path_str):
            m = MagicMock()
            m.exists.return_value = 'genome2' in str(path_str)  # genome1 never becomes ready
            return m

        processed = []

        def fake_run_subprocess(cmd):
            processed.append(_target_genome_of(cmd))
            return self._success_proc()

        with patch('bioinformatics_tools.workflow_tools.workflow.Path', side_effect=fake_path), \
             patch.object(wf, '_start_background_subprocess', return_value=fake_background), \
             patch.object(wf, '_run_subprocess', side_effect=fake_run_subprocess):
            wf._run_pipeline_batch_sequential('margie_sb', smk_config, genome_files,
                                              cache_map_fn=lambda s: {}, mode='dev')

        assert processed == ['genome2']
        assert wf.report.status.indicates_success


# ---------------------------------------------------------------------------
# do_margie_sb
# ---------------------------------------------------------------------------

class TestDoMargieSb:

    def test_genome_cache_map_covers_every_phase4_8_tool(self, wf):
        """do_margie_sb should build a real per-tool cache map (not the old
        2-file whole-genome marker) and pass it through as cache_map_fn --
        every phase4-8 tool present, with the known multi-file tools
        (tigrfam, phobius, envelope, interpro) carrying their full file set."""
        wf.conf.get = MagicMock(side_effect=lambda key, default=None: {
            'input': '/genomes/genome1.fa',
            'main_database': '/tmp/test-margie.db',
        }.get(key, default))

        captured = {}

        def fake_run_pipeline_batch_sequential(key_name, smk_config, genome_files, **kwargs):
            captured['cache_map_fn'] = kwargs['cache_map_fn']

        with patch('bioinformatics_tools.workflow_tools.workflow.discover_genomes',
                   return_value={'genome1': '/genomes/genome1.fa'}), \
             patch.object(wf, '_run_pipeline_batch_sequential', side_effect=fake_run_pipeline_batch_sequential):
            wf.do_margie_sb(mode='dev')

        cache_map = captured['cache_map_fn']('genome1')

        expected_tools = {
            'cog', 'pfam', 'merops', 'tcdb', 'uniprot', 'kegg', 'eggnog', 'dbcan',
            'pgap', 'geneprop', 'operon', 'tmbed', 'signalp6', 'deepsig', 'psortb',
            'signalp4', 'tigrfam', 'phobius', 'envelope', 'interpro',
        }
        assert set(cache_map.keys()) == expected_tools
        # quast/gtdbtk/rasttk deliberately excluded -- see _genome_cache_map's docstring.
        assert 'quast' not in cache_map
        assert 'gtdbtk' not in cache_map
        assert 'rasttk' not in cache_map

        assert len(cache_map['cog']) == 2  # results.tsv + token
        assert len(cache_map['tigrfam']) == 3  # results.tsv + domtbl.out + token
        assert len(cache_map['phobius']) == 3  # results.tsv + top1.tsv + token
        assert len(cache_map['envelope']) == 3  # results.tsv + summary.tsv + token
        assert len(cache_map['interpro']) == 38  # (results.tsv + token) + 18 dbs * 2 files each

        for paths in cache_map.values():
            for path in paths:
                assert path.endswith(('.tsv', '.tkn', '.out'))


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
