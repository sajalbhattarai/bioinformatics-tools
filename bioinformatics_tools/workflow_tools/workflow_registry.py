"""
Central registry of all available workflows.

Defines all Snakemake workflows that can be executed via the dane_wf CLI.
Each workflow is registered as a WorkflowKey with metadata for execution,
frontend display, and configuration.
"""
import logging
import os

from bioinformatics_tools.workflow_tools.models import WorkflowKey
from bioinformatics_tools.workflow_tools.workflow_helpers import WORKFLOW_PATH_DEFAULTS

LOGGER = logging.getLogger(__name__)

# Root for the per-user writable stores a run appends to (OCC operon reference,
# fingerprint databases, genome pool, scoring archive). {user} is filled in by
# resolve_user_paths(); read-only reference data stays shared.
MARGIE_DEPOT_ROOT = os.environ.get('MARGIE_SHARED_ROOT', '').rstrip('/')
MARGIE_USER_ROOT = f'{MARGIE_DEPOT_ROOT}/users/{{user}}'


def resolve_user_paths(params: list[dict], username: str | None) -> list[dict]:
    """Returns *params* with {user} in path defaults filled in for *username*.

    *username* is the cluster username. Without one, the shared legacy root is
    used (and logged) rather than creating a literal "{user}" directory.
    """
    if not params:
        return params
    if username:
        target = f'{MARGIE_DEPOT_ROOT}/users/{username}'
    else:
        LOGGER.warning('No cluster username available to resolve per-user store paths; '
                       'falling back to the shared depot root')
        target = MARGIE_DEPOT_ROOT
    token = f'{MARGIE_DEPOT_ROOT}/users/{{user}}'
    out = []
    for p in params:
        d = p.get('default')
        if isinstance(d, str) and token in d:
            p = {**p, 'default': d.replace(token, target)}
        out.append(p)
    return out


# System-wide required parameters for cluster execution
# These are needed for ANY workflow running via SLURM, not workflow-specific
REQUIRED_SYSTEM_PARAMS = [
    {
        'param': 'compute.cluster_default.account',
        'default': None,
        'description': 'SLURM account for job submission (REQUIRED for cluster execution)',
        'type': 'string',
        'required': True
    },
    {
        'param': 'compute.cluster_default.partition',
        'default': 'cpu',
        'description': 'SLURM partition/queue for job submission',
        'type': 'string'
    },
    {
        'param': 'compute.cluster_default.default_runtime',
        'default': 240,
        'description': 'Default runtime limit in minutes for SLURM jobs',
        'type': 'int'
    },
    {
        'param': 'compute.cluster_default.default_mem_mb',
        'default': 4000,
        'description': 'Default memory limit in MB for SLURM jobs',
        'type': 'int'
    },
    {
        'param': 'compute.cluster_default.max_jobs',
        'default': 5,
        'description': 'Maximum number of concurrent SLURM jobs',
        'type': 'int'
    },
    {
        # Walltime of the driver job that runs Snakemake and queues every rule
        # job (default_runtime caps those). Adjustable so it can be lowered to
        # fit before a maintenance reservation, which blocks longer requests.
        'param': 'compute.cluster_default.driver_walltime',
        'default': '3-00:00:00',
        'description': ('SLURM walltime for the workflow driver job (D-HH:MM:SS or HH:MM:SS). '
                        'Must exceed the whole run. Lower it to fit before a maintenance '
                        'outage, otherwise the job stays PENDING until the outage ends.'),
        'type': 'string'
    },
]


def workflow_path_params(wf_id: str, include_sif: bool = True, include_db_root: bool = True,
                         supports_batch_input: bool = True) -> list[dict]:
    """Returns the per-workflow root-path settings (sif_path, db_root, input_path,
    output_path), namespaced as e.g. 'margie_sb.sif_path'.

    sif_path/db_root are included only for workflows with local SIF lookup
    (local_sif_only) or a db_root fallback (supports_db_root).
    """
    defaults = WORKFLOW_PATH_DEFAULTS.get(wf_id, {})
    params = []
    if include_sif:
        params.append({
            'param': f'{wf_id}.sif_path',
            'default': defaults.get('sif_path', '~/.cache/bioinformatics-tools'),
            'description': f'Root folder containing all Apptainer .sif container files for {wf_id} (flat, e.g. sif_path/cog.sif)',
            'type': 'path'
        })
    if include_db_root:
        params.append({
            'param': f'{wf_id}.db_root',
            'default': defaults.get('db_root', ''),
            'description': f'Root folder containing per-tool database subfolders for {wf_id} (e.g. db_root/cog/) — not every tool needs one',
            'type': 'path'
        })
    input_desc = (f'Default genome file or folder of genomes for {wf_id}, used when a run does not specify one explicitly'
                  if supports_batch_input else
                  f'Default single genome file for {wf_id}, used when a run does not specify one explicitly')
    params.append({
        'param': f'{wf_id}.input_path',
        'default': defaults.get('input_path', ''),
        'description': input_desc,
        'type': 'path'
    })
    params.append({
        'param': f'{wf_id}.output_path',
        'default': defaults.get('output_path', ''),
        'description': f'Default output root for {wf_id}, used when a run does not specify one explicitly',
        'type': 'path'
    })
    return params


# 'database': the tool reads a reference database from db.<key>, else
# <margie_sb.db_root>/<key> (margie_sb.smk's db_path calls).
MARGIE_SB_PHASED_TOOLS = [
    {'key': 'quast', 'label': 'QUAST', 'phase': 1, 'sif': 'quast.sif', 'purpose': 'Assembly quality metrics'},
    {'key': 'gtdbtk', 'label': 'GTDB-Tk', 'phase': 2, 'sif': 'gtdbtk.sif', 'purpose': 'Taxonomic assignment', 'database': True},
    {'key': 'rasttk', 'label': 'RASTtk', 'phase': 3, 'sif': 'rasttk.sif', 'purpose': 'Core annotation stage gate', 'database': True},
    {'key': 'cog', 'label': 'COG', 'phase': 4, 'sif': 'cog.sif', 'purpose': 'Functional category annotation', 'database': True},
    {'key': 'kegg', 'label': 'KEGG', 'phase': 4, 'sif': 'kegg.sif', 'purpose': 'Pathway annotation', 'database': True},
    {'key': 'eggnog', 'label': 'eggNOG', 'phase': 4, 'sif': 'eggnog.sif', 'purpose': 'Orthology annotation', 'database': True},
    {'key': 'uniprot', 'label': 'UniProt', 'phase': 4, 'sif': 'uniprot.sif', 'purpose': 'Protein annotation', 'database': True},
    {'key': 'pfam', 'label': 'Pfam', 'phase': 4, 'sif': 'pfam.sif', 'purpose': 'Protein family annotation', 'database': True},
    {'key': 'tigrfam', 'label': 'TIGRFAM', 'phase': 4, 'sif': 'tigrfam.sif', 'purpose': 'Protein family annotation', 'database': True},
    {'key': 'merops', 'label': 'MEROPS', 'phase': 4, 'sif': 'merops.sif', 'purpose': 'Protease annotation', 'database': True},
    {'key': 'tcdb', 'label': 'TCDB', 'phase': 4, 'sif': 'tcdb.sif', 'purpose': 'Transporter annotation', 'database': True},
    {'key': 'dbcan', 'label': 'dbCAN', 'phase': 4, 'sif': 'dbcan.sif', 'purpose': 'CAZyme annotation', 'database': True},
    {'key': 'pgap', 'label': 'PGAP', 'phase': 4, 'sif': 'pgap.sif', 'purpose': 'Genome annotation support', 'database': True},
    {'key': 'geneprop', 'label': 'GeneProp', 'phase': 4, 'sif': 'geneprop.sif', 'purpose': 'Gene property annotation (after TIGRFAM)', 'database': True},
    {'key': 'interpro', 'label': 'InterPro', 'phase': 4, 'sif': 'interpro.sif', 'purpose': 'Domain/signature annotation', 'database': True},
    {'key': 'operon', 'label': 'Operon', 'phase': 5, 'sif': 'operon.sif', 'purpose': 'Operon prediction'},
    {'key': 'phobius', 'label': 'Phobius', 'phase': 6, 'sif': 'phobius.sif', 'purpose': 'Signal peptide and topology prediction'},
    {'key': 'tmbed', 'label': 'TMbed', 'phase': 6, 'sif': 'tmbed.sif', 'purpose': 'Transmembrane topology prediction', 'database': True},
    # SignalP 4 and 6 run from RCAC's own modules (envmodules:, not container:
    # in margie_sb.smk), so neither has a .sif to check or set up.
    {'key': 'signalp6', 'label': 'SignalP6', 'phase': 6, 'sif': 'signalp6.sif', 'purpose': 'Signal peptide prediction (HPC module, no envelope dependency)', 'uses_container': False},
    {'key': 'envelope', 'label': 'Envelope', 'phase': 7, 'sif': 'envelope.sif', 'purpose': 'Gram envelope type inference'},
    {'key': 'psortb', 'label': 'PSORTb', 'phase': 8, 'sif': 'psortb.sif', 'purpose': 'Subcellular localization'},
    {'key': 'deepsig', 'label': 'DeepSig', 'phase': 8, 'sif': 'deepsig.sif', 'purpose': 'Signal peptide prediction'},
    {'key': 'signalp4', 'label': 'SignalP4', 'phase': 8, 'sif': 'signalP4.sif', 'purpose': 'Signal peptide prediction (HPC module; kept until SignalP6 is confirmed)', 'uses_container': False},
    {'key': 'consolidation', 'label': 'Consolidation', 'phase': 9, 'sif': 'consolidation.sif', 'purpose': 'Consolidate upstream outputs', 'uses_container': False},
    {'key': 'labeling', 'label': 'Labeling', 'phase': 10, 'sif': 'labeling.sif', 'purpose': 'Label assignment', 'uses_container': False},
    # Host-side Python like consolidation/labeling (margie_sb.smk's run_fingerprint and
    # update_fingerprint_database have no container:), so no .sif is checked or listed.
    {'key': 'fingerprint', 'label': 'Fingerprint', 'phase': 11, 'sif': 'fingerprint.sif', 'purpose': 'Feature fingerprinting', 'uses_container': False},
    {'key': 'scoring_heuristic', 'label': 'Scoring', 'phase': 12, 'sif': 'scoring-heuristic.sif', 'purpose': 'Confidence scoring (statistical model)', 'uses_container': False},
    {'key': 'fingerprint_database', 'label': 'Fingerprint Database', 'phase': 13, 'sif': 'fingerprint-database.sif', 'purpose': 'Fingerprint DB stage', 'uses_container': False},
    {'key': 'ani', 'label': 'ANI', 'phase': 14, 'sif': 'ani.sif', 'purpose': 'Average nucleotide identity'},
    {'key': 'aai', 'label': 'AAI', 'phase': 14, 'sif': 'aai.sif', 'purpose': 'Average amino acid identity'},
    {'key': 'closest', 'label': 'Closest', 'phase': 14, 'sif': 'closest.sif', 'purpose': 'Closest genome matching'},
    # No rule in margie_sb.smk runs these yet, so no .sif is checked or listed.
    {'key': 'mauve', 'label': 'Mauve', 'phase': 14, 'sif': 'mauve.sif', 'purpose': 'Whole-genome synteny/collinear blocks vs. closest organisms', 'uses_container': False},
    {'key': 'synteny', 'label': 'Synteny', 'phase': 14, 'sif': 'synteny.sif', 'purpose': 'Synteny calculation', 'uses_container': False},
    {'key': 'evidence', 'label': 'Evidence', 'phase': 15, 'sif': 'evidence.sif', 'purpose': 'Per-gene annotation evidence reports (build-gene-report.py)', 'uses_container': False},
    {'key': 'llm', 'label': 'LLM', 'phase': 15, 'sif': 'llm.sif', 'purpose': 'LLM confidence scoring layer (GPU, ROCm/MI210)'},
]


def _margie_sb_default_threads(tool_key: str) -> int:
    if tool_key == 'gtdbtk':
        return 64
    if tool_key in {'kegg', 'eggnog'}:
        return 16
    return 8


def _margie_sb_default_mem_mb(tool_key: str) -> int:
    if tool_key == 'gtdbtk':
        return 460000
    if tool_key == 'interpro':
        return 48000
    if tool_key == 'eggnog':
        return 64000
    if tool_key == 'kegg':
        return 16000
    if tool_key == 'dbcan':
        return 16000
    if tool_key == 'pgap':
        return 12000
    if tool_key == 'tmbed':
        return 32000
    return 4000


def _margie_sb_default_runtime(tool_key: str) -> int:
    if tool_key in {'interpro'}:
        return 300
    return 240


def margie_sb_sif_files(selected_tool_keys: set[str] | None = None) -> list[tuple]:
    """Returns the SIF entries to validate for a margie_sb run.

    Skips tools with uses_container=False and, when selected_tool_keys is given,
    tools outside that subset.
    """
    return [
        (tool['sif'], 'latest') for tool in MARGIE_SB_PHASED_TOOLS
        if tool.get('uses_container', True)
        and (selected_tool_keys is None or tool['key'] in selected_tool_keys)
    ]


def _margie_sb_location_params() -> list[dict]:
    """Returns the params for each tool's database location (when not
    <db_root>/<tool>) and margie-build's location on the cluster."""
    params = [{
        'param': f"db.{tool['key']}",
        'default': '',
        'description': f"{tool['label']} database folder; empty means <Databases folder>/{tool['key']}",
        'type': 'path',
    } for tool in MARGIE_SB_PHASED_TOOLS if tool.get('database')]
    params.append({
        'param': 'margie_sb.build_repo',
        'default': '',
        'description': "margie-build's folder on the cluster, which builds containers and databases; "
                       "empty means the lab's copy on depot, else ~/margie-build",
        'type': 'path',
    })
    return params


def _margie_sb_tool_params() -> list[dict]:
    params: list[dict] = []
    for tool in MARGIE_SB_PHASED_TOOLS:
        key = tool['key']
        label = tool['label']
        phase = tool['phase']
        params.extend([
            {
                'param': f'margie_sb.{key}.threads',
                'default': _margie_sb_default_threads(key),
                'description': f'Phase {phase}: thread count for {label}',
                'type': 'int'
            },
            {
                'param': f'margie_sb.{key}.partition',
                'default': 'highmem' if key == 'gtdbtk' else 'cpu',
                'description': f'Phase {phase}: SLURM partition for {label}',
                'type': 'string'
            },
            {
                'param': f'margie_sb.{key}.mem_mb',
                'default': _margie_sb_default_mem_mb(key),
                'description': f'Phase {phase}: memory limit (MB) for {label}',
                'type': 'int'
            },
            {
                'param': f'margie_sb.{key}.runtime',
                'default': _margie_sb_default_runtime(key),
                'description': f'Phase {phase}: runtime limit (minutes) for {label}',
                'type': 'int'
            },
            {
                'param': f'margie_sb.{key}.db',
                'default': '',
                'description': f'Phase {phase}: database path for {label}',
                'type': 'path'
            },
            {
                'param': f'margie_sb.{key}.sif',
                'default': tool['sif'],
                'description': f'Phase {phase}: SIF filename for {label} under sif_path',
                'type': 'string'
            },
        ])
    return params


WORKFLOWS: dict[str, WorkflowKey] = {
    'example': WorkflowKey(
        cmd_identifier='example',
        snakemake_file='example.smk',
        other=[''],
        sif_files=[
            ('prodigal.sif', '2.6.3-v1.0'),
        ],
        label='Example',
        description='Simple test workflow for development',
        full_description='A minimal workflow for testing the pipeline infrastructure.',
    ),
    'margie': WorkflowKey(
        cmd_identifier='margie',
        snakemake_file='margie.smk',
        other=[''],
        sif_files=[
            ('prodigal.sif', '2.6.3-v1.0'),
            ('pfam_scan_light', 'latest'),
            ('cogclassifier', 'latest'),
            ('kofam_scan_light_bsp', 'latest'),
            ('opr-dev', '1.13'),
            ('run_dbcan_light', '4.2.0')
        ],
        label='Margie',
        description='Full annotation pipeline (Prodigal, Pfam, COG)',
        full_description='Comprehensive microbial genome annotation workflow that combines gene prediction with functional annotation. Runs Prodigal for open reading frame prediction, Pfam for protein family identification, and COGclassifier for functional categorization. Results are automatically loaded into a SQLite database for downstream analysis.',
        tools=[
            {
                'name': 'Prodigal',
                'purpose': 'Gene prediction and ORF identification',
                'version': '2.6.3',
                'output': 'GFF3 file with predicted genes and protein sequences (FAA)'
            },
            {
                'name': 'Pfam_scan',
                'purpose': 'Protein family and domain annotation',
                'version': 'latest',
                'output': 'CSV file with Pfam domain hits'
            },
            {
                'name': 'COGclassifier',
                'purpose': 'Functional categorization using COG database',
                'version': 'latest',
                'output': 'TSV files with COG classifications and category counts'
            },
            {
                'name': 'Kegg',
                'purpose': 'TODO',
                'version': 'latest',
                'output': 'kegg.tkn # TODO'
            },
            {
                'name': 'UniOp',
                'purpose': 'TODO',
                'version': 'latest',
                'output': '# TODO'
            },
            {
                'name': 'run_dbCAN',
                'purpose': 'TODO',
                'version': 'latest',
                'output': '# TODO'
            }
        ],
        configurable_params=[
            # Prodigal configuration (rule run_prodigal)
            {
                'param': 'prodigal.threads',
                'default': 1,
                'description': 'Number of threads for Prodigal',
                'type': 'int'
            },
            {
                'param': 'prodigal.mem_mb',
                'default': 2048,
                'description': 'Memory limit in MB for Prodigal',
                'type': 'int'
            },
            {
                'param': 'prodigal.runtime',
                'default': 30,
                'description': 'Runtime limit in minutes for Prodigal',
                'type': 'int'
            },
            # Pfam configuration (rule run_pfam)
            {
                'param': 'pfam.threads',
                'default': 4,
                'description': 'Number of threads for Pfam scan',
                'type': 'int'
            },
            {
                'param': 'pfam.mem_mb',
                'default': 4000,
                'description': 'Memory limit in MB for Pfam scan',
                'type': 'int'
            },
            {
                'param': 'pfam.runtime',
                'default': 240,
                'description': 'Runtime limit in minutes for Pfam scan',
                'type': 'int'
            },
            {
                'param': 'pfam.db',
                'default': '/depot/lindems/data/Databases/pfam',
                'description': 'Path to Pfam-A HMM database',
                'type': 'path'
            },
            # UniOp configuration (rule run_uniop)
            {
                'param': 'uniop.output_dir',
                'default': 'uniop',
                'description': 'Output directory for uniop - requires a directory, not file',
                'type': 'str'
            },
            # COG configuration (rule run_cog)
            {
                'param': 'cog.threads',
                'default': 4,
                'description': 'Number of threads for COGclassifier',
                'type': 'int'
            },
            {
                'param': 'cog.mem_mb',
                'default': 8192,
                'description': 'Memory limit in MB for COGclassifier',
                'type': 'int'
            },
            {
                'param': 'cog.runtime',
                'default': 120,
                'description': 'Runtime limit in minutes for COGclassifier',
                'type': 'int'
            },
            {
                'param': 'cog.db',
                'default': '/depot/lindems/data/Databases/cog/',
                'description': 'Path to COG database directory',
                'type': 'path'
            },
            {
                'param': 'cog.outdir',
                'default': 'cog',
                'description': 'Output directory for COG results',
                'type': 'string'
            },
            # dbCAN configuration (rule run_dbcan)
            {
                'param': 'dbcan.threads',
                'default': 4,
                'description': 'Number of threads for dbCAN',
                'type': 'int'
            },
            {
                'param': 'dbcan.mem_mb',
                'default': 7984,
                'description': 'Memory limit in MB for dbCAN',
                'type': 'int'
            },
            {
                'param': 'dbcan.runtime',
                'default': 180,
                'description': 'Runtime limit in minutes for dbCAN',
                'type': 'int'
            },
            {
                'param': 'dbcan.db',
                'default': '/depot/lindems/data/Databases/cazyme/db',
                'description': 'Path to dbCAN database directory',
                'type': 'path'
            },
            {
                'param': 'dbcan.output_dir',
                'default': 'dbcan',
                'description': 'Output directory for dbCAN results',
                'type': 'string'
            }
        ],
        database_deps=[
            'Pfam-A HMM profiles',
            'COG functional database',
            'SQLite results database'
        ],
        docs_url=None
    ),
    'margie_sb': WorkflowKey(
        cmd_identifier='margie_sb',
        snakemake_file='margie_sb.smk',
        other=[''],
        sif_files=margie_sb_sif_files(),
        local_sif_only=True,
        supports_batch_input=True,
        supports_db_root=True,
        label='MARGIE (SB)',
        description='Custom phased MARGIE workflow by sajalbhattarai',
        full_description='Phased MARGIE(SB) workflow wiring with full container inventory and per-tool resource/path configuration. Phase1 and phase2 can continue with warnings if a tool fails, phase3+ enforce strict dependency gates before downstream phases proceed.',
        tools=[
            {
                'key': tool['key'],
                'phase': tool['phase'],
                'name': tool['label'],
                'purpose': f"Phase {tool['phase']}: {tool['purpose']}",
                'version': 'latest',
                'output': 'Phase report in margie_sb/phases plus stage-specific outputs in later functional commits'
            }
            for tool in MARGIE_SB_PHASED_TOOLS
        ],
        configurable_params=[
            {
                'param': 'margie_sb.default_threads',
                'default': 8,
                'description': 'Default thread count for MARGIE(SB) tools unless overridden per tool',
                'type': 'int'
            },
            {
                'param': 'margie_sb.default_mem_mb',
                'default': 4000,
                'description': 'Default memory limit in MB for MARGIE(SB) tools unless overridden per tool',
                'type': 'int'
            },
            {
                'param': 'margie_sb.default_runtime',
                'default': 240,
                'description': 'Default runtime limit in minutes for MARGIE(SB) tools unless overridden per tool',
                'type': 'int'
            },
            # Folder for each user's growing databases and run archives, copied
            # from the bases on the first run and tracked by api/services/user_stores.py.
            # Empty: /scratch/<cluster>/<cluster user>/margie-2026.
            {
                'param': 'margie_sb.stores_root',
                'default': '',
                'description': 'Working folder for your databases and run archives (empty: margie-2026 in your scratch)',
                'type': 'path'
            },
            {
                'param': 'margie_sb.backup_root',
                'default': f'{MARGIE_DEPOT_ROOT}/databases/margie-generated-databases' if MARGIE_DEPOT_ROOT else '',
                'description': 'Where "Back up to depot" puts copies of your databases, beside the bases every user starts from',
                'type': 'path'
            },
            {
                'param': 'margie_sb.phase4.max_parallel_tools',
                'default': 4,
                'description': 'Max parallel phase4 annotation tools per genome',
                'type': 'int'
            },
            {
                'param': 'margie_sb.phase8.max_parallel_tools',
                'default': 4,
                'description': 'Max parallel phase8 (envelope-dependent localization) tools per genome',
                'type': 'int'
            },
            {
                'param': 'margie_sb.phase1.partition',
                'default': 'cpu',
                'description': 'Default SLURM partition for phase1 tools',
                'type': 'string'
            },
            {
                'param': 'margie_sb.phase2.partition',
                'default': 'cpu',
                'description': 'Default SLURM partition for phase2 tools',
                'type': 'string'
            },
            {
                'param': 'margie_sb.phase3.partition',
                'default': 'cpu',
                'description': 'Default SLURM partition for phase3 tools',
                'type': 'string'
            },
            {
                'param': 'margie_sb.phase3.max_parallel_genomes',
                'default': 1,
                'description': 'Max parallel phase3 RASTtk genomes (default 1 to protect BV-BRC service)',
                'type': 'int'
            },
            {
                'param': 'margie_sb.phase4.partition',
                'default': 'cpu',
                'description': 'Default SLURM partition for phase4 tools',
                'type': 'string'
            },
            {
                'param': 'margie_sb.phase5.partition',
                'default': 'cpu',
                'description': 'Default SLURM partition for phase5 tools',
                'type': 'string'
            },
            {
                'param': 'margie_sb.phase6.partition',
                'default': 'cpu',
                'description': 'Default SLURM partition for phase6 tools',
                'type': 'string'
            },
            {
                'param': 'margie_sb.phase7.partition',
                'default': 'cpu',
                'description': 'Default SLURM partition for phase7 tools',
                'type': 'string'
            },
            {
                'param': 'margie_sb.phase8.partition',
                'default': 'cpu',
                'description': 'Default SLURM partition for phase8 tools',
                'type': 'string'
            },
            {
                'param': 'margie_sb.phase9.partition',
                'default': 'cpu',
                'description': 'Default SLURM partition for phase9 tools',
                'type': 'string'
            },
            {
                'param': 'margie_sb.phase10.partition',
                'default': 'cpu',
                'description': 'Default SLURM partition for phase10 tools',
                'type': 'string'
            },
            {
                'param': 'margie_sb.phase11.partition',
                'default': 'cpu',
                'description': 'Default SLURM partition for phase11 tools',
                'type': 'string'
            },
            {
                'param': 'margie_sb.phase12.partition',
                'default': 'cpu',
                'description': 'Default SLURM partition for phase12 tools',
                'type': 'string'
            },
            {
                'param': 'margie_sb.phase13.partition',
                'default': 'cpu',
                'description': 'Default SLURM partition for phase13 tools',
                'type': 'string'
            },
            {
                'param': 'margie_sb.phase14.partition',
                'default': 'cpu',
                'description': 'Default SLURM partition for phase14 tools',
                'type': 'string'
            },
            {
                'param': 'margie_sb.phase15.partition',
                'default': 'cpu',
                'description': 'Default SLURM partition for phase15 tools',
                'type': 'string'
            },
        ] + _margie_sb_tool_params() + _margie_sb_location_params(),
        database_deps=[
            'Input FASTA file',
            'Configurable db_root with per-tool db overrides',
            'Configurable sif_path with per-tool sif filename overrides'
        ],
        docs_url=None
    ),
    'selftest': WorkflowKey(
        cmd_identifier='selftest',
        snakemake_file='selftest.smk',
        other=[''],
        sif_files=[],
        label='Self Test',
        description='Quick validation test (no containers)',
        full_description='Lightweight test workflow that validates SSH, Snakemake, and database caching without using containers. Useful for verifying the pipeline infrastructure is working correctly.',
    ),
}


def get_workflow(name: str) -> WorkflowKey | None:
    """
    Get a workflow definition by name.

    Args:
        name: The workflow identifier (e.g., 'margie', 'example')

    Returns:
        WorkflowKey if found, None otherwise
    """
    return WORKFLOWS.get(name)


def list_workflows() -> list[str]:
    """
    List all registered workflow names.

    Returns:
        List of workflow identifiers
    """
    return list(WORKFLOWS.keys())
