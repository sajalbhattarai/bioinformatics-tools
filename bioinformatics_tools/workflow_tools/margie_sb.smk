"""
MARGIE(SB) phased workflow scaffold.

This workflow declares all requested container phases and validates configurable
SIF/DB roots. It currently emits phase reports and markers while preserving
strict phase order and warning/continue behavior where requested.
"""
import json
import os
import sys
from pathlib import Path

# Add current directory to path to import workflow_helpers
sys.path.insert(0, os.path.dirname(workflow.snakefile))
from workflow_helpers import rc, fixed_path, db_token, db_path, sif_dir


INPUT_FASTA = rc('input_fasta', config=config)
MAIN_DATABASE = rc('main_database', config=config)
SIF_DIR = sif_dir(config=config)

MARGIE_SB_READY = fixed_path('margie_sb/margie_sb_ready.txt', config=config)
MARGIE_SB_TOKEN = db_token('margie_sb', config=config)

PHASE1_DONE = fixed_path('margie_sb/phases/phase1.done', config=config)
PHASE2_DONE = fixed_path('margie_sb/phases/phase2.done', config=config)
PHASE3_DONE = fixed_path('margie_sb/phases/phase3.done', config=config)
PHASE4_DONE = fixed_path('margie_sb/phases/phase4.done', config=config)
PHASE5_DONE = fixed_path('margie_sb/phases/phase5.done', config=config)
PHASE6_DONE = fixed_path('margie_sb/phases/phase6.done', config=config)
PHASE7_DONE = fixed_path('margie_sb/phases/phase7.done', config=config)
PHASE8_DONE = fixed_path('margie_sb/phases/phase8.done', config=config)
PHASE9_DONE = fixed_path('margie_sb/phases/phase9.done', config=config)
PHASE10_DONE = fixed_path('margie_sb/phases/phase10.done', config=config)
PHASE11_DONE = fixed_path('margie_sb/phases/phase11.done', config=config)
PHASE12_DONE = fixed_path('margie_sb/phases/phase12.done', config=config)
PHASE13_DONE = fixed_path('margie_sb/phases/phase13.done', config=config)


PHASES = {
    'phase1': {
        'description': 'Quality checks before classification',
        'allow_continue_on_failure': True,
        'tools': [
            {'key': 'quast', 'sif': 'quast.sif'},
            {'key': 'checkm', 'sif': 'checkm.sif'},
        ],
    },
    'phase2': {
        'description': 'Organism classification (bacteria/archaea)',
        'allow_continue_on_failure': True,
        'tools': [
            {'key': 'gtdbtk', 'sif': 'gtdbtk.sif'},
            {'key': 'classifier', 'sif': 'classifier.sif'},
        ],
    },
    'phase3': {
        'description': 'Core annotation gate (must pass before phase4)',
        'allow_continue_on_failure': False,
        'tools': [
            {'key': 'rasttk', 'sif': 'rasttk.sif'},
        ],
    },
    'phase4': {
        'description': 'Annotation tools (parallelizable by config)',
        'allow_continue_on_failure': False,
        'tools': [
            {'key': 'cog', 'sif': 'cog.sif'},
            {'key': 'kegg', 'sif': 'kegg.sif'},
            {'key': 'eggnog', 'sif': 'eggnog.sif'},
            {'key': 'uniprot', 'sif': 'uniprot.sif'},
            {'key': 'pfam', 'sif': 'pfam.sif'},
            {'key': 'tigrfam', 'sif': 'tigrfam.sif'},
            {'key': 'merops', 'sif': 'merops.sif'},
            {'key': 'tcdb', 'sif': 'tcdb.sif'},
            {'key': 'dbcan', 'sif': 'dbcan.sif'},
            {'key': 'pgap', 'sif': 'pgap.sif'},
            {'key': 'geneprop', 'sif': 'geneprop.sif'},
            {'key': 'interpro', 'sif': 'interpro.sif'},
            {'key': 'operon', 'sif': 'operon.sif'},
            {'key': 'tmbed', 'sif': 'tmbed.sif'},
        ],
    },
    'phase5': {
        'description': 'Gram envelope type inference',
        'allow_continue_on_failure': False,
        'tools': [
            {'key': 'envelope', 'sif': 'envelope.sif'},
        ],
    },
    'phase6': {
        'description': 'Localization and topology',
        'allow_continue_on_failure': False,
        'tools': [
            {'key': 'tmhmm', 'sif': 'tmhmm.sif'},
            {'key': 'tatfinder', 'sif': 'tatfinder.sif'},
            {'key': 'phobius', 'sif': 'phobius.sif'},
            {'key': 'psortb', 'sif': 'psortb.sif'},
            {'key': 'deepsig', 'sif': 'deepsig.sif'},
            {'key': 'signalp4', 'sif': 'signalP4.sif'},
        ],
    },
    'phase7': {
        'description': 'Consolidates upstream outputs',
        'allow_continue_on_failure': False,
        'tools': [
            {'key': 'consolidation', 'sif': 'consolidation.sif'},
        ],
    },
    'phase8': {
        'description': 'Labeling',
        'allow_continue_on_failure': False,
        'tools': [
            {'key': 'labeling', 'sif': 'labeling.sif'},
        ],
    },
    'phase9': {
        'description': 'Fingerprint generation',
        'allow_continue_on_failure': False,
        'tools': [
            {'key': 'fingerprint', 'sif': 'fingerprint.sif'},
        ],
    },
    'phase10': {
        'description': 'Scoring heuristic',
        'allow_continue_on_failure': False,
        'tools': [
            {'key': 'scoring_heuristic', 'sif': 'scoring-heuristic.sif'},
        ],
    },
    'phase11': {
        'description': 'Fingerprint database',
        'allow_continue_on_failure': False,
        'tools': [
            {'key': 'fingerprint_database', 'sif': 'fingerprint-database.sif'},
        ],
    },
    'phase12': {
        'description': 'Synteny calculations',
        'allow_continue_on_failure': False,
        'tools': [
            {'key': 'ani', 'sif': 'ani.sif'},
            {'key': 'aai', 'sif': 'aai.sif'},
            {'key': 'closest', 'sif': 'closest.sif'},
            {'key': 'synteny', 'sif': 'synteny.sif'},
        ],
    },
    'phase13': {
        'description': 'LLM-based analysis',
        'allow_continue_on_failure': False,
        'tools': [
            {'key': 'llm', 'sif': 'llm.sif'},
        ],
    },
}


def _tool_value(tool_key: str, field: str, default):
    return rc(
        f'margie_sb.{tool_key}.{field}',
        rc(f'margie_sb.default_{field}', default, config=config),
        config=config,
    )


def _phase_partition(phase_key: str) -> str:
    return rc(
        f'margie_sb.{phase_key}.partition',
        rc('compute.cluster_default.partition', 'cpu', config=config),
        config=config,
    )


def _validate_phase(phase_key: str, output_path: str) -> None:
    phase = PHASES[phase_key]
    allow_continue = phase['allow_continue_on_failure']

    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    warnings: list[str] = []
    missing_required: list[str] = []
    report: list[dict] = []

    for tool in phase['tools']:
        tool_key = tool['key']
        sif_name = _tool_value(tool_key, 'sif', tool['sif'])
        resolved_sif = Path(SIF_DIR) / sif_name
        db_value = _tool_value(tool_key, 'db', db_path(tool_key, config=config))

        exists = resolved_sif.exists()
        report.append(
            {
                'phase': phase_key,
                'tool': tool_key,
                'sif': sif_name,
                'resolved_sif': str(resolved_sif),
                'sif_exists': exists,
                'db': db_value,
                'threads': _tool_value(tool_key, 'threads', 8),
                'mem_mb': _tool_value(tool_key, 'mem_mb', 4000),
                'runtime': _tool_value(tool_key, 'runtime', 120),
                'partition': _tool_value(tool_key, 'partition', _phase_partition(phase_key)),
            }
        )

        if not exists:
            warning = f'[WARN] Missing container for {tool_key}: {resolved_sif}'
            warnings.append(warning)
            if not allow_continue:
                missing_required.append(warning)

    with out.open('w', encoding='utf-8') as handle:
        handle.write(f'phase={phase_key}\n')
        handle.write(f"description={phase['description']}\n")
        handle.write(f'allow_continue_on_failure={allow_continue}\n')
        for warning in warnings:
            handle.write(f'{warning}\n')
        handle.write(json.dumps(report, indent=2))
        handle.write('\n')

    if missing_required:
        raise ValueError('\n'.join(missing_required))


rule all:
    input:
        MARGIE_SB_TOKEN


rule margie_sb_scaffold:
    """Validate base inputs and write workflow-level readiness metadata."""
    input:
        INPUT_FASTA
    output:
        MARGIE_SB_READY
    threads: rc('margie_sb.default_threads', 8, config=config)
    resources:
        mem_mb=rc('margie_sb.default_mem_mb', 4000, config=config),
        runtime=rc('margie_sb.default_runtime', 120, config=config),
        slurm_partition=rc('compute.cluster_default.partition', 'cpu', config=config)
    run:
        out = Path(output[0])
        out.parent.mkdir(parents=True, exist_ok=True)
        scaffold = {
            'input_fasta': str(input[0]),
            'main_database': MAIN_DATABASE,
            'sif_root': SIF_DIR,
            'db_root': rc('db_root', '/depot/lindems/data/Databases', config=config),
            'phase4_max_parallel_tools': rc('margie_sb.phase4.max_parallel_tools', 4, config=config),
            'phase6_max_parallel_tools': rc('margie_sb.phase6.max_parallel_tools', 4, config=config),
            'phase_dependency_note': 'geneprop should run after tigrfam in phase4 execution logic',
        }
        with out.open('w', encoding='utf-8') as handle:
            handle.write(json.dumps(scaffold, indent=2))
            handle.write('\n')


rule margie_sb_phase1:
    input:
        MARGIE_SB_READY
    output:
        PHASE1_DONE
    run:
        _validate_phase('phase1', output[0])


rule margie_sb_phase2:
    input:
        PHASE1_DONE
    output:
        PHASE2_DONE
    run:
        _validate_phase('phase2', output[0])


rule margie_sb_phase3:
    input:
        PHASE2_DONE
    output:
        PHASE3_DONE
    run:
        _validate_phase('phase3', output[0])


rule margie_sb_phase4:
    input:
        PHASE3_DONE
    output:
        PHASE4_DONE
    run:
        _validate_phase('phase4', output[0])


rule margie_sb_phase5:
    input:
        PHASE4_DONE
    output:
        PHASE5_DONE
    run:
        _validate_phase('phase5', output[0])


rule margie_sb_phase6:
    input:
        PHASE5_DONE
    output:
        PHASE6_DONE
    run:
        _validate_phase('phase6', output[0])


rule margie_sb_phase7:
    input:
        PHASE6_DONE
    output:
        PHASE7_DONE
    run:
        _validate_phase('phase7', output[0])


rule margie_sb_phase8:
    input:
        PHASE7_DONE
    output:
        PHASE8_DONE
    run:
        _validate_phase('phase8', output[0])


rule margie_sb_phase9:
    input:
        PHASE8_DONE
    output:
        PHASE9_DONE
    run:
        _validate_phase('phase9', output[0])


rule margie_sb_phase10:
    input:
        PHASE9_DONE
    output:
        PHASE10_DONE
    run:
        _validate_phase('phase10', output[0])


rule margie_sb_phase11:
    input:
        PHASE10_DONE
    output:
        PHASE11_DONE
    run:
        _validate_phase('phase11', output[0])


rule margie_sb_phase12:
    input:
        PHASE11_DONE
    output:
        PHASE12_DONE
    run:
        _validate_phase('phase12', output[0])


rule margie_sb_phase13:
    input:
        PHASE12_DONE
    output:
        PHASE13_DONE
    run:
        _validate_phase('phase13', output[0])


rule margie_sb_finalize:
    """Create workflow token for cache compatibility in later iterations."""
    input:
        PHASE13_DONE
    output:
        MARGIE_SB_TOKEN
    shell:
        """
        mkdir -p $(dirname {output})
        touch {output}
        """
