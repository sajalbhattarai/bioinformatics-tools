"""
MARGIE(SB) scaffold workflow.

This is the initial registration target for the custom Snakemake pipeline.
It is intentionally minimal and will be expanded in later commits.
"""
import os
import sys

# Add current directory to path to import workflow_helpers
sys.path.insert(0, os.path.dirname(workflow.snakefile))
from workflow_helpers import rc, fixed_path, db_token, db_path, sif_dir


INPUT_FASTA = rc('input_fasta', config=config)
MAIN_DATABASE = rc('main_database', config=config)
SIF_DIR = sif_dir(config=config)
KEGG_DB = db_path('kegg', config=config)

# Placeholder outputs for first-commit scaffold.
MARGIE_SB_READY = fixed_path('margie_sb/margie_sb_ready.txt', config=config)
MARGIE_SB_TOKEN = db_token('margie_sb', config=config)


rule all:
    input:
        MARGIE_SB_TOKEN


rule margie_sb_scaffold:
    """Scaffold rule that validates config wiring and creates a readiness marker."""
    input:
        INPUT_FASTA
    output:
        MARGIE_SB_READY
    threads: rc('margie_sb.threads', 1, config=config)
    resources:
        mem_mb=rc('margie_sb.mem_mb', 1024, config=config),
        runtime=rc('margie_sb.runtime', 5, config=config)
    shell:
        """
        mkdir -p $(dirname {output})
        printf "MARGIE(SB) scaffold ready\ninput=%s\nmain_database=%s\nsif_dir=%s\nkegg_db=%s\n" \
            "{input}" "{MAIN_DATABASE}" "{SIF_DIR}" "{KEGG_DB}" > {output}
        """


rule margie_sb_finalize:
    """Create workflow token for cache compatibility in later iterations."""
    input:
        MARGIE_SB_READY
    output:
        MARGIE_SB_TOKEN
    shell:
        """
        mkdir -p $(dirname {output})
        touch {output}
        """
