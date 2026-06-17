"""
Central registry of all available workflows.

Defines all Snakemake workflows that can be executed via the dane_wf CLI.
Each workflow is registered as a WorkflowKey with metadata for execution,
frontend display, and configuration.
"""
from bioinformatics_tools.workflow_tools.models import WorkflowKey


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
        'default': 30,
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
]


MARGIE_SB_PHASED_TOOLS = [
    {'key': 'quast', 'label': 'QUAST', 'phase': 1, 'sif': 'quast.sif', 'purpose': 'Assembly quality metrics'},
    {'key': 'checkm', 'label': 'CheckM', 'phase': 1, 'sif': 'checkm.sif', 'purpose': 'Genome quality/completeness checks'},
    {'key': 'gtdbtk', 'label': 'GTDB-Tk', 'phase': 2, 'sif': 'gtdbtk.sif', 'purpose': 'Taxonomic assignment'},
    {'key': 'classifier', 'label': 'Classifier', 'phase': 2, 'sif': 'classifier.sif', 'purpose': 'Bacteria/archaea classification'},
    {'key': 'rasttk', 'label': 'RASTtk', 'phase': 3, 'sif': 'rasttk.sif', 'purpose': 'Core annotation stage gate'},
    {'key': 'cog', 'label': 'COG', 'phase': 4, 'sif': 'cog.sif', 'purpose': 'Functional category annotation'},
    {'key': 'kegg', 'label': 'KEGG', 'phase': 4, 'sif': 'kegg.sif', 'purpose': 'Pathway annotation'},
    {'key': 'eggnog', 'label': 'eggNOG', 'phase': 4, 'sif': 'eggnog.sif', 'purpose': 'Orthology annotation'},
    {'key': 'uniprot', 'label': 'UniProt', 'phase': 4, 'sif': 'uniprot.sif', 'purpose': 'Protein annotation'},
    {'key': 'pfam', 'label': 'Pfam', 'phase': 4, 'sif': 'pfam.sif', 'purpose': 'Protein family annotation'},
    {'key': 'tigrfam', 'label': 'TIGRFAM', 'phase': 4, 'sif': 'tigrfam.sif', 'purpose': 'Protein family annotation'},
    {'key': 'merops', 'label': 'MEROPS', 'phase': 4, 'sif': 'merops.sif', 'purpose': 'Protease annotation'},
    {'key': 'tcdb', 'label': 'TCDB', 'phase': 4, 'sif': 'tcdb.sif', 'purpose': 'Transporter annotation'},
    {'key': 'dbcan', 'label': 'dbCAN', 'phase': 4, 'sif': 'dbcan.sif', 'purpose': 'CAZyme annotation'},
    {'key': 'pgap', 'label': 'PGAP', 'phase': 4, 'sif': 'pgap.sif', 'purpose': 'Genome annotation support'},
    {'key': 'geneprop', 'label': 'GeneProp', 'phase': 4, 'sif': 'geneprop.sif', 'purpose': 'Gene property annotation (after TIGRFAM)'},
    {'key': 'interpro', 'label': 'InterPro', 'phase': 4, 'sif': 'interpro.sif', 'purpose': 'Domain/signature annotation'},
    {'key': 'operon', 'label': 'Operon', 'phase': 4, 'sif': 'operon.sif', 'purpose': 'Operon prediction'},
    {'key': 'tmbed', 'label': 'TMbed', 'phase': 4, 'sif': 'tmbed.sif', 'purpose': 'Membrane topology support'},
    {'key': 'envelope', 'label': 'Envelope', 'phase': 5, 'sif': 'envelope.sif', 'purpose': 'Gram envelope type inference'},
    {'key': 'tmhmm', 'label': 'TMHMM', 'phase': 6, 'sif': 'tmhmm.sif', 'purpose': 'Transmembrane helix prediction'},
    {'key': 'tatfinder', 'label': 'TatFinder', 'phase': 6, 'sif': 'tatfinder.sif', 'purpose': 'Tat signal detection'},
    {'key': 'phobius', 'label': 'Phobius', 'phase': 6, 'sif': 'phobius.sif', 'purpose': 'Signal peptide and topology prediction'},
    {'key': 'psortb', 'label': 'PSORTb', 'phase': 6, 'sif': 'psortb.sif', 'purpose': 'Subcellular localization'},
    {'key': 'deepsig', 'label': 'DeepSig', 'phase': 6, 'sif': 'deepsig.sif', 'purpose': 'Signal peptide prediction'},
    {'key': 'signalp4', 'label': 'SignalP4', 'phase': 6, 'sif': 'signalP4.sif', 'purpose': 'Signal peptide prediction'},
    {'key': 'consolidation', 'label': 'Consolidation', 'phase': 7, 'sif': 'consolidation.sif', 'purpose': 'Consolidate upstream outputs'},
    {'key': 'labeling', 'label': 'Labeling', 'phase': 8, 'sif': 'labeling.sif', 'purpose': 'Label assignment'},
    {'key': 'fingerprint', 'label': 'Fingerprint', 'phase': 9, 'sif': 'fingerprint.sif', 'purpose': 'Feature fingerprinting'},
    {'key': 'scoring_heuristic', 'label': 'Scoring Heuristic', 'phase': 10, 'sif': 'scoring-heuristic.sif', 'purpose': 'Heuristic scoring'},
    {'key': 'fingerprint_database', 'label': 'Fingerprint Database', 'phase': 11, 'sif': 'fingerprint-database.sif', 'purpose': 'Fingerprint DB stage'},
    {'key': 'ani', 'label': 'ANI', 'phase': 12, 'sif': 'ani.sif', 'purpose': 'Average nucleotide identity'},
    {'key': 'aai', 'label': 'AAI', 'phase': 12, 'sif': 'aai.sif', 'purpose': 'Average amino acid identity'},
    {'key': 'closest', 'label': 'Closest', 'phase': 12, 'sif': 'closest.sif', 'purpose': 'Closest genome matching'},
    {'key': 'synteny', 'label': 'Synteny', 'phase': 12, 'sif': 'synteny.sif', 'purpose': 'Synteny calculation'},
    {'key': 'llm', 'label': 'LLM', 'phase': 13, 'sif': 'llm.sif', 'purpose': 'LLM-based analysis'},
]


def _margie_sb_default_threads(tool_key: str) -> int:
    if tool_key in {'kegg', 'eggnog'}:
        return 16
    return 8


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
                'default': 'cpu',
                'description': f'Phase {phase}: SLURM partition for {label}',
                'type': 'string'
            },
            {
                'param': f'margie_sb.{key}.mem_mb',
                'default': 4000,
                'description': f'Phase {phase}: memory limit (MB) for {label}',
                'type': 'int'
            },
            {
                'param': f'margie_sb.{key}.runtime',
                'default': 120,
                'description': f'Phase {phase}: runtime limit (minutes) for {label}',
                'type': 'int'
            },
            {
                'param': f'margie_sb.{key}.db',
                'default': f'/depot/lindems/data/Databases/{key}',
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
        sif_files=[(tool['sif'], 'latest') for tool in MARGIE_SB_PHASED_TOOLS],
        label='MARGIE (SB)',
        description='Custom phased MARGIE workflow by sajalbhattarai',
        full_description='Phased MARGIE(SB) workflow wiring with full container inventory and per-tool resource/path configuration. Phase1 and phase2 can continue with warnings if a tool fails, phase3+ enforce strict dependency gates before downstream phases proceed.',
        tools=[
            {
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
                'default': 120,
                'description': 'Default runtime limit in minutes for MARGIE(SB) tools unless overridden per tool',
                'type': 'int'
            },
            {
                'param': 'db_root',
                'default': '/depot/lindems/data/Databases',
                'description': 'Root directory for per-tool databases (all MARGIE(SB) stages)',
                'type': 'path'
            },
            {
                'param': 'sif_path',
                'default': '~/.cache/bioinformatics-tools',
                'description': 'Directory containing Apptainer SIF files for all MARGIE(SB) stages',
                'type': 'path'
            },
            {
                'param': 'margie_sb.phase4.max_parallel_tools',
                'default': 4,
                'description': 'Max parallel phase4 annotation tools per genome',
                'type': 'int'
            },
            {
                'param': 'margie_sb.phase6.max_parallel_tools',
                'default': 4,
                'description': 'Max parallel phase6 localization tools per genome',
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
        ] + _margie_sb_tool_params(),
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
