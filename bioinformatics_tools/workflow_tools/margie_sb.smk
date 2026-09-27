# MARGIE single-bacterium (margie_sb) Snakemake workflow: per-genome QC, taxonomy,
# gene calling, functional annotation, localization and database loading.

import os
import re
import sys

# Makes workflow_helpers importable from this directory.
sys.path.insert(0, os.path.dirname(workflow.snakefile))
from workflow_helpers import rc, rc_bool, fixed_path, sif_path, db_path, db_token, discover_genomes, genome_calls, default_store_root, get_workflow_prefix_for, get_container_outputs_prefix_for
from load_to_db import PIPELINE_VERSION

WORKFLOW_DIR = os.path.dirname(workflow.snakefile)
LOAD_SCRIPT = os.path.join(WORKFLOW_DIR, "load_to_db.py")
_repo_root = os.path.dirname(os.path.dirname(WORKFLOW_DIR))
_default_python = os.path.join(_repo_root, ".venv", "bin", "python")
LOADER_PYTHON = os.environ.get("MARGIE_PYTHON", _default_python)
if not os.path.exists(LOADER_PYTHON):
    raise ValueError(
        f"Required Python interpreter not found: {LOADER_PYTHON}. "
        "Use the repo-local .venv (uv sync) or set MARGIE_PYTHON to a valid path."
    )
# Appends the genome-level ENVELOPE_* columns to a phase8 tool's results.tsv.
ENRICH_SCRIPT = os.path.join(WORKFLOW_DIR, "enrich_with_envelope.py")
# Host-side processing scripts for signalp6/signalp4, which run as HPC
# environment modules and have no container entrypoint of their own.
SIGNALP6_SCRIPT = os.path.join(WORKFLOW_DIR, "process_signalp6.py")
SIGNALP4_SCRIPT = os.path.join(WORKFLOW_DIR, "process_signalp4.py")


# Store root when the config names none: the user's scratch, never the depot bases.
STORE_ROOT = rc('margie_sb.stores_root', '', config=config) or default_store_root()
# The group's base databases (MARGIE_SHARED_ROOT on the cluster); without one, the user's own stores.
_SHARED = os.environ.get('MARGIE_SHARED_ROOT', '').rstrip('/')
BASES = f'{_SHARED}/databases/margie-generated-databases' if _SHARED else STORE_ROOT


def _resolve_cfg_path(preferred_key: str, legacy_key: str, default: str) -> str:
    """Returns a config path from the namespaced key, falling back to the legacy key."""
    value = rc(preferred_key, rc(legacy_key, default, config=config), config=config)
    return str(value).strip()


def _resolve_shared_dir(preferred_key: str, legacy_key: str, default: str) -> str:
    """Returns a directory path from config without a trailing slash."""
    return _resolve_cfg_path(preferred_key, legacy_key, default).rstrip("/")


def _resolve_shared_file(preferred_key: str, legacy_key: str, default: str, canonical_name: str) -> str:
    """Returns a file path from config, appending canonical_name when only a directory is given."""
    raw = _resolve_cfg_path(preferred_key, legacy_key, default)
    trimmed = raw.rstrip("/")
    leaf = os.path.basename(trimmed)
    if raw.endswith("/") or not leaf or "." not in leaf:
        return f"{trimmed}/{canonical_name}"
    return raw

# ---- path definitions (single source of truth for all rule paths) ----

# input_fasta is one FASTA or a directory of them; discover_genomes() maps
# genome stem -> file. GENOME_PREFIX carries the literal {genome} wildcard,
# so every rule runs once per discovered genome.
MAIN_DATABASE = rc('main_database', config=config)

GENOMES = discover_genomes(rc('input_fasta', config=config))
if not GENOMES:
    raise ValueError(f"No genome files found for input_fasta={rc('input_fasta', config=config)!r}")

GENOME_PREFIX = get_workflow_prefix_for('{genome}', config=config)
# Scratch root that each tool's container writes its raw output to; output_dir
# (GENOME_PREFIX) holds only the tracked results.tsv files and db tokens.
CONTAINER_OUTPUTS_PREFIX = get_container_outputs_prefix_for('{genome}', config=config)

_OUTPUT_ROOT = rc('output_dir', '', config=config).rstrip('/')
MIN_SLURM_RUNTIME_MINUTES = 240


def runtime_min(key: str, default: int, config=None) -> int:
    """Returns the configured rule runtime in minutes, clamped to at least 4 hours."""
    minutes = rc(key, default, config=config)
    try:
        minutes = int(minutes)
    except (TypeError, ValueError):
        minutes = int(default)
    return max(MIN_SLURM_RUNTIME_MINUTES, minutes)

# Account-wide mutex directory serializing run_rasttk's BV-BRC submissions.
# It lives in the store root, not a per-run path (which can inherit a stale lock from a
# staged older run) nor $HOME (not fully bind-mounted into the rasttk container).
RASTTK_BVBRC_LOCK = _resolve_shared_dir(
    'rasttk.bvbrc_lock_dir', 'rasttk_bvbrc_lock_dir',
    f'{STORE_ROOT}/rasttk_bvbrc.lock',
)
os.makedirs(os.path.dirname(RASTTK_BVBRC_LOCK), exist_ok=True)

# Quast outputs
QUAST_RESULTS = f"{GENOME_PREFIX}quast/quast.tsv"
QUAST_TOKEN = f"{GENOME_PREFIX}quast/quast_db.tkn"

# Batch staging and aggregation paths for QUAST.
QUAST_BATCH_PREFIX = f"{_OUTPUT_ROOT}/original_container_outputs/quast" if _OUTPUT_ROOT else "original_container_outputs/quast"
QUAST_BATCH_STAGE_DIR = f"{QUAST_BATCH_PREFIX}/stage"
QUAST_BATCH_OUTPUT_DIR = f"{QUAST_BATCH_PREFIX}/container_outputs"
QUAST_BATCH_DONE = f"{QUAST_BATCH_PREFIX}/quast_batch.done"

# GTDB-Tk outputs
GTDBTK_RESULTS = f"{GENOME_PREFIX}gtdbtk/gtdbtk_results.tsv"
# Genetic code for phase3 (RASTtk); not loaded into the database.
GTDBTK_TRANSLATION_TABLE = f"{GENOME_PREFIX}gtdbtk/translation_table.tsv"
GTDBTK_TOKEN = f"{GENOME_PREFIX}gtdbtk/gtdbtk_db.tkn"
GTDBTK_COMPUTE_TOKEN = f"{GENOME_PREFIX}gtdbtk/gtdbtk_compute.tkn"

# Batch paths for GTDB-Tk: one classify_wf over all genomes (shared DB warm-up),
# then split back into per-genome GTDBTK_RESULTS / GTDBTK_TRANSLATION_TABLE.
GTDBTK_BATCH_PREFIX = f"{_OUTPUT_ROOT}/original_container_outputs/gtdbtk" if _OUTPUT_ROOT else "original_container_outputs/gtdbtk"
GTDBTK_BATCH_STAGE_DIR = f"{GTDBTK_BATCH_PREFIX}/stage"
GTDBTK_BATCH_OUTPUT_DIR = f"{GTDBTK_BATCH_PREFIX}/container_outputs"
GTDBTK_BATCH_RESULTS = f"{GTDBTK_BATCH_PREFIX}/gtdbtk_results.tsv"
GTDBTK_BATCH_TRANSLATION_TABLE = f"{GTDBTK_BATCH_PREFIX}/gtdbtk.translation_table_summary.tsv"
GTDBTK_BATCH_DONE = f"{GTDBTK_BATCH_PREFIX}/gtdbtk_batch.done"

# RASTtk outputs. rast.faa/.gff are inputs to phase4 and operon, so they stay
# in output_dir; the other gene_calls files are flattened into rasttk/ untracked.
RASTTK_RESULTS = f"{GENOME_PREFIX}rasttk/rast.tsv"
RASTTK_TOKEN = f"{GENOME_PREFIX}rasttk/rasttk_db.tkn"
RASTTK_COMPUTE_TOKEN = f"{GENOME_PREFIX}rasttk/rasttk_compute.tkn"
RASTTK_FAA = f"{GENOME_PREFIX}rasttk/rast.faa"
RASTTK_GFF = f"{GENOME_PREFIX}rasttk/rast.gff"

# Per-protein cache (protein_cache.py): workflow.py writes <tool>/protein-cache/
# novel.faa with the uncached proteins; the rule annotates only those and
# PC_MERGE (a POSIX sh script) merges its rows with the cached ones. Without
# novel.faa the rule reads rast.faa and the merge is a plain copy.
import protein_cache as _protein_cache
PC_MERGE = _protein_cache.install_merge_script(_OUTPUT_ROOT or os.path.join(WORKFLOW_DIR, '.protein-cache-local'))


def pc_dir(tool):
    # Returns the tool's protein-cache directory for a genome.
    return lambda wildcards: f"{GENOME_PREFIX}{tool}/{_protein_cache.PC_DIR}".format(genome=wildcards.genome)


def pc_faa(tool):
    # Returns the tool's novel.faa when present, otherwise rast.faa.
    def _faa(wildcards):
        novel = f"{GENOME_PREFIX}{tool}/{_protein_cache.PC_DIR}/novel.faa".format(genome=wildcards.genome)
        return novel if os.path.exists(novel) else RASTTK_FAA.format(genome=wildcards.genome)
    return _faa

# Each genome's domain, genetic code and gene caller (genome_calls). With
# GTDB-Tk off, both come from margie_sb.genome_info; a genome missing either
# is called by Prodigal, which writes RASTtk's layout so later rules are shared.
RUN_GTDBTK = rc_bool('run_gtdbtk', True, config=config)
GENOME_CALLS = genome_calls(GENOMES, config)
RASTTK_GENOMES = sorted(g for g, c in GENOME_CALLS.items() if c['gene_caller'] == 'rasttk')
PRODIGAL_GENOMES = sorted(g for g, c in GENOME_CALLS.items() if c['gene_caller'] == 'prodigal')


def _one_of(names):
    """Returns a {genome} wildcard regex matching exactly these genomes (nothing if empty)."""
    return '|'.join(re.escape(n) for n in names) if names else '(?!)'


# Per-genome domain (GTDBTK_domain) and genetic code (translation_table), using
# GTDB-Tk's column names; written from GTDB-Tk or genome_info, config wins.
GENOME_INFO = f"{GENOME_PREFIX}genome_info/genome_info.tsv"

# Prodigal's own container outputs; the gene calls land in RASTTK_* above.
PRODIGAL_CONTAINER_OUTPUTS = f"{CONTAINER_OUTPUTS_PREFIX}prodigal"

# Phase4: functional annotation (12 tools). Each entrypoint takes
# -i <faa> -o <root> -d <db> -t <threads> --organism-name {genome} --domain <d>
# and writes <root>/{genome}/processed/<tool>_results.tsv. mem_mb defaults
# track each tool's database size. margie_sb_phase4_slot caps concurrent
# phase4 tools (set from margie_sb.phase4.max_parallel_tools by workflow.py).
PHASE4_TOOLS = [
    "pgap", "tigrfam", "uniprot", "pfam", "kegg", "eggnog", "cog",
    "merops", "tcdb", "dbcan", "geneprop", "interpro",
]
PHASE4_RESULTS = {t: f"{GENOME_PREFIX}{t}/{t}_results.tsv" for t in PHASE4_TOOLS}
PHASE4_TOKENS = {t: f"{GENOME_PREFIX}{t}/{t}_db.tkn" for t in PHASE4_TOOLS}
PHASE4_COMPUTE_TOKENS = {t: f"{GENOME_PREFIX}{t}/{t}_compute.tkn" for t in PHASE4_TOOLS}
# Raw tigrfam hmmscan domtblout, an input to run_geneprop (--tigrfam-domtbl).
TIGRFAM_DOMTBL = f"{GENOME_PREFIX}tigrfam/tigrfam_domtbl.out"

# InterPro per-database outputs. PHASE4_RESULTS['interpro'] is the unified
# table (member database in INTERPRO_analysis); the split TSVs are loaded too.
# All member analyses available in this InterProScan install (display name ->
# basename); Phobius/SignalP/TMHMM are excluded as they run standalone.
INTERPRO_ALL_ANALYSES = {
    "AntiFam": "antifam", "CDD": "cdd", "Coils": "coils", "FunFam": "funfam",
    "Gene3D": "gene3d", "Hamap": "hamap", "MobiDBLite": "mobidb", "NCBIfam": "ncbifam",
    "PANTHER": "panther", "Pfam": "pfam", "PIRSF": "pirsf", "PIRSR": "pirsr",
    "PRINTS": "prints", "ProSitePatterns": "prosite_patterns", "ProSiteProfiles": "prosite_profiles",
    "SFLD": "sfld", "SMART": "smart", "SUPERFAMILY": "superfamily",
}

# Active subset from interpro.analyses; the default four are prokaryote-focused,
# avoid overlap with standalone Pfam/TIGRFAM, and keep run_interpro's memory
# under the per-node limit. The per-db paths and load rule derive from it.
_INTERPRO_ACTIVE_NAMES = rc('interpro.analyses', ["Hamap", "NCBIfam", "CDD", "PIRSF"], config=config)
_unknown = [n for n in _INTERPRO_ACTIVE_NAMES if n not in INTERPRO_ALL_ANALYSES]
if _unknown:
    raise ValueError(
        f"interpro.analyses: unknown analysis name(s) {_unknown}; "
        f"must be a subset of {sorted(INTERPRO_ALL_ANALYSES)}"
    )
INTERPRO_ANALYSIS_TO_BASENAME = {name: INTERPRO_ALL_ANALYSES[name] for name in _INTERPRO_ACTIVE_NAMES}
INTERPRO_DEFAULT_APPS = ",".join(INTERPRO_ANALYSIS_TO_BASENAME.keys())
INTERPRO_DB_BASENAMES = list(INTERPRO_ANALYSIS_TO_BASENAME.values())
INTERPRO_PERDB_RESULTS = {
    db: f"{GENOME_PREFIX}interpro/interpro_{db}_results.tsv" for db in INTERPRO_DB_BASENAMES
}
INTERPRO_PERDB_TOKEN_PATTERN = f"{GENOME_PREFIX}interpro/interpro_{{db}}_db.tkn"

# Phase5: operon prediction (UniOP) from rast.faa plus rast.gff (-i/-g); no database.
OPERON_RESULTS = f"{GENOME_PREFIX}operon/operon_results.tsv"
OPERON_TOKEN = f"{GENOME_PREFIX}operon/operon_db.tkn"
OPERON_COMPUTE_TOKEN = f"{GENOME_PREFIX}operon/operon_compute.tkn"

# Phase6: envelope-independent localization/topology (phobius, tmbed, signalp6),
# each reading a single FAA.
PHOBIUS_RESULTS = f"{GENOME_PREFIX}phobius/phobius_results.tsv"
PHOBIUS_TOKEN = f"{GENOME_PREFIX}phobius/phobius_db.tkn"
PHOBIUS_COMPUTE_TOKEN = f"{GENOME_PREFIX}phobius/phobius_compute.tkn"
# Per-protein Phobius summary; not loaded into the database.
PHOBIUS_TOP1 = f"{GENOME_PREFIX}phobius/phobius_top1.tsv"

# TMbed: deep-learning transmembrane predictor using the ProtT5 model in db/tmbed; runs on CPU.
TMBED_RESULTS = f"{GENOME_PREFIX}tmbed/tmbed_results.tsv"
TMBED_TOKEN = f"{GENOME_PREFIX}tmbed/tmbed_db.tkn"
TMBED_COMPUTE_TOKEN = f"{GENOME_PREFIX}tmbed/tmbed_compute.tkn"

# SignalP 6.0 runs as an HPC environment module (envmodules:, not container:).
# --format none avoids per-protein plot files named after the whole FASTA header.
# The processing script is overridable via signalp6.process_script.
SIGNALP6_RESULTS = f"{GENOME_PREFIX}signalp6/signalp6_results.tsv"
SIGNALP6_TOKEN = f"{GENOME_PREFIX}signalp6/signalp6_db.tkn"
SIGNALP6_COMPUTE_TOKEN = f"{GENOME_PREFIX}signalp6/signalp6_compute.tkn"
SIGNALP6_PROCESS_SCRIPT = rc('signalp6.process_script', SIGNALP6_SCRIPT, config=config)

# Phase7: envelope type inference (monoderm vs diderm). -i is the genome's
# output root, which the entrypoint searches for the tigrfam/pgap/pfam/uniprot
# results; those four are declared as inputs only to order the DAG.
ENVELOPE_RESULTS = f"{GENOME_PREFIX}envelope/envelope_results.tsv"
ENVELOPE_TOKEN = f"{GENOME_PREFIX}envelope/envelope_db.tkn"
ENVELOPE_COMPUTE_TOKEN = f"{GENOME_PREFIX}envelope/envelope_compute.tkn"
# Genome-level envelope decision (always one row), read by the phase8 tools.
ENVELOPE_SUMMARY = f"{GENOME_PREFIX}envelope/envelope_summary.tsv"

# Phase8: envelope-dependent localization (deepsig, psortb, signalp4).
# DeepSig's -k maps ENVELOPE_SUMMARY's envelope_type to GRAM-/GRAM+/ARCH.
DEEPSIG_RESULTS = f"{GENOME_PREFIX}deepsig/deepsig_results.tsv"
DEEPSIG_TOKEN = f"{GENOME_PREFIX}deepsig/deepsig_db.tkn"
DEEPSIG_COMPUTE_TOKEN = f"{GENOME_PREFIX}deepsig/deepsig_compute.tkn"

# PSORTb v3: like deepsig, but -k takes n|p|a; the entrypoint tolerates
# PSORTb's non-zero exits on warnings.
PSORTB_RESULTS = f"{GENOME_PREFIX}psortb/psortb_results.tsv"
PSORTB_TOKEN = f"{GENOME_PREFIX}psortb/psortb_db.tkn"
PSORTB_COMPUTE_TOKEN = f"{GENOME_PREFIX}psortb/psortb_compute.tkn"

# SignalP 4.1 runs as an HPC environment module (command `signalp`). Its -t has
# no archaea option, so archaea map to gram- in run_signalp4. The processing
# script is overridable via signalp4.process_script.
SIGNALP4_RESULTS = f"{GENOME_PREFIX}signalp4/signalp4_results.tsv"
SIGNALP4_TOKEN = f"{GENOME_PREFIX}signalp4/signalp4_db.tkn"
SIGNALP4_COMPUTE_TOKEN = f"{GENOME_PREFIX}signalp4/signalp4_compute.tkn"
SIGNALP4_PROCESS_SCRIPT = rc('signalp4.process_script', SIGNALP4_SCRIPT, config=config)

# Phase9 (consolidation): host-side scripts in consolidation/ (detect-columns,
# merge-all-columns, filter-no-stat); no container.
CONSOLIDATION_SCRIPTS_DIR = os.path.join(WORKFLOW_DIR, "consolidation")
CONSOLIDATION_DETECTED_COLUMNS = f"{GENOME_PREFIX}consolidation/detected-columns.json"
CONSOLIDATION_MERGED = f"{GENOME_PREFIX}consolidation/consolidated-merged-all-columns.tsv"
CONSOLIDATION_MANIFEST = f"{GENOME_PREFIX}consolidation/manifest.tsv"
CONSOLIDATION_NO_STAT = f"{GENOME_PREFIX}consolidation/consolidated-no-stat.tsv"
# COMPUTE_TOKEN marks files written (read by later phases); TOKEN marks
# computed and loaded into main_database (requested by rule all).
CONSOLIDATION_COMPUTE_TOKEN = f"{GENOME_PREFIX}consolidation/consolidation_compute.tkn"
CONSOLIDATION_TOKEN = f"{GENOME_PREFIX}consolidation/consolidation_db.tkn"

# Phase10 (labeling): scripts in labeling/. assign-canonical-label.py reads the
# full merged table (it needs InterPro and score columns); the other three each
# read the labeled output plus the merged table.
LABELING_SCRIPTS_DIR = os.path.join(WORKFLOW_DIR, "labeling")
LABELING_LABELED = f"{GENOME_PREFIX}labeling/labeled-genes.tsv"
LABELING_EC_CONSENSUS = f"{GENOME_PREFIX}labeling/labeled-genes-ec-consensus.tsv"
LABELING_OPERON_INFO = f"{GENOME_PREFIX}labeling/labeled-genes-operon-info.tsv"
# Collapses the TIGRFAM/PGAP/NCBIfam HMM cluster into one C1 slot and adds
# COG/KEGG corroboration signals.
LABELING_CLUSTER_AGREEMENT = f"{GENOME_PREFIX}labeling/labeled-genes-cluster-agreement.tsv"
LABELING_COMPUTE_TOKEN = f"{GENOME_PREFIX}labeling/labeling_compute.tkn"
LABELING_TOKEN = f"{GENOME_PREFIX}labeling/labeling_db.tkn"

# Phase11 (scoring): host-side scripts in scoring/ computing the hierarchy and
# confidence tiers, the C1-C4 components and the blended confidence_score.
SCORING_SCRIPTS_DIR = os.path.join(WORKFLOW_DIR, "scoring")
# Persistent cross-run OCC operon reference (user's copy on scratch).
C3_REFERENCE_PKL = _resolve_shared_file(
    'margie_sb.operon_database.occ_reference_pkl',
    'operon_database.occ_reference_pkl',
    f'{STORE_ROOT}/operon-database/occ_reference.pkl',
    'occ_reference.pkl',
)
SCORING_HIERARCHY_TIER = f"{GENOME_PREFIX}scoring/scored-labeled-genes-annotation-tool-tier.tsv"
SCORING_CONFIDENCE_TIER = f"{GENOME_PREFIX}scoring/scored-labeled-genes-annotation-ec-tier.tsv"
SCORING_C1 = f"{GENOME_PREFIX}scoring/scored-labeled-genes-c1-tool-coverage.tsv"
SCORING_C2 = f"{GENOME_PREFIX}scoring/scored-labeled-genes-c2-operon-probability.tsv"
SCORING_C3 = f"{GENOME_PREFIX}scoring/scored-labeled-genes-c3-operonic-context-confidence.tsv"
SCORING_C4 = f"{GENOME_PREFIX}scoring/scored-labeled-genes-c4-ec-agreement.tsv"
SCORING_CONFIDENCE_FINAL = f"{GENOME_PREFIX}scoring/scored-labeled-genes-confidence-final.tsv"
SCORING_FINAL_ANNOTATION_WITH_CONFIDENCE = f"{GENOME_PREFIX}scoring/FINAL_ANNOTATION_WITH_CONFIDENCE.tsv"
SCORING_OCC_REFERENCE_TOKEN = f"{GENOME_PREFIX}scoring/occ_reference_updated.tkn"
SCORING_COMPUTE_TOKEN = f"{GENOME_PREFIX}scoring/scoring_compute.tkn"
SCORING_TOKEN = f"{GENOME_PREFIX}scoring/scoring_db.tkn"

# Archive of every run's final scoring table. C3 depends on a growing OCC
# reference, so the database keeps only the latest scores and each run's
# tables are snapshotted here for history.
SCORING_HISTORICAL_PATH = _resolve_shared_dir(
    'margie_sb.scoring_results_historical.path',
    'scoring_results_historical.path',
    f'{STORE_ROOT}/scoring-archive',
)
# One archive folder per run, named after the run's timestamped output directory.
_RUN_TIMESTAMP = os.path.basename(_OUTPUT_ROOT) if _OUTPUT_ROOT else 'adhoc'
SCORING_ARCHIVE_DIR = f"{SCORING_HISTORICAL_PATH}/{_RUN_TIMESTAMP}"
SCORING_ARCHIVE_TOKEN = f"{GENOME_PREFIX}scoring/scoring_archived.tkn"

# Per-organism export: <store>/final-tables/<organism>/FINAL_ANNOTATION_WITH_CONFIDENCE.tsv
FINAL_TABLES_DEPOT_PATH = _resolve_shared_dir(
    'margie_sb.final_tables_depot.path',
    'final_tables_depot.path',
    f'{STORE_ROOT}/final-tables',
)
FINAL_TABLES_DEPOT_TOKEN = f"{GENOME_PREFIX}scoring/final_tables_depot.tkn"

# ---- post-scoring report figures ----
# Figures and TSVs written to the run's output tree from finished scoring outputs
# and the read-only operon reference; the rules never fail, so they cannot block
# scoring. Per-organism: <genome>/scoring/figures/; pangenome: <run>/scoring/figures/global/.
REPORT_FIGURES_SCRIPTS_DIR = os.path.join(SCORING_SCRIPTS_DIR, "analysis", "report_figures")
REPORT_FIGURES_OPERON_DB = _resolve_shared_file(
    'margie_sb.report_figures.operon_db',
    'report_figures.operon_db',
    # Read-only, so the depot base is a safe fallback.
    f'{BASES}/fingerprint-database/operon-fingerprint-database-label-ordered.tsv',
    'operon-fingerprint-database-label-ordered.tsv',
)
REPORT_FIGURES_ORGANISM_DIR = f"{GENOME_PREFIX}scoring/figures"
REPORT_FIGURES_ORGANISM_TOKEN = f"{GENOME_PREFIX}scoring/report_figures.tkn"
REPORT_FIGURES_GLOBAL_DIR = f"{_OUTPUT_ROOT}/scoring/figures/global"
REPORT_FIGURES_GLOBAL_TOKEN = f"{_OUTPUT_ROOT}/scoring/figures/report_figures_global.tkn"

# Interactive genome/operon viewer: one self-contained HTML file per organism,
# written at the organism top level (reorganize_outputs.py keeps GENOME_VIEWER_NAME).
VIZ_SCRIPTS_DIR = os.path.join(WORKFLOW_DIR, "viz")
GENOME_VIEWER_NAME = "FINAL_GENOME_VIEWER.html"
GENOME_VIEWER_HTML = f"{GENOME_PREFIX}{GENOME_VIEWER_NAME}"
GENOME_VIEWER_CIRCULAR_PNG = f"{GENOME_PREFIX}scoring/figures/{{genome}}_circular.png"
GENOME_VIEWER_TOKEN = f"{GENOME_PREFIX}scoring/genome_viewer.tkn"
# Optional full operon atlas (every multi-gene operon) under
# <genome>/scoring/figures/complete-organism-operon-diagrams/; off by default.
COMPLETE_OPERON_MAP_TOKEN = f"{GENOME_PREFIX}scoring/figures/complete_operon_map.tkn"


def _home_config_flag(section, key):
    """Returns a boolean setting from ~/.config/bioinformatics-tools/config.yaml.

    Returns False when the file or key is missing or unreadable.
    """
    try:
        import yaml
        p = os.path.join(os.path.expanduser("~"), ".config",
                         "bioinformatics-tools", "config.yaml")
        with open(p) as fh:
            cfg = yaml.safe_load(fh) or {}
        val = (cfg.get(section) or {}).get(key, False)
        if isinstance(val, str):
            return val.strip().lower() not in ("false", "0", "no", "off", "")
        return bool(val)
    except Exception:
        return False


# Default for run_full_operon_map from the home config; the run config overrides it.
_FULL_OPERON_MAP_DEFAULT = _home_config_flag("margie_sb", "run_full_operon_map")

# Target of a non-blocking SLURM job that snapshots margie.db by pipeline version.
SQLITE_SNAPSHOT_ROOT = _resolve_shared_dir(
    'margie_sb.sqlite_pipeline_snapshot.path',
    'sqlite_pipeline_snapshot.path',
    f'{STORE_ROOT}/sqlite/snapshots',
)
SQLITE_SNAPSHOT_VERSION_DIR = f"{SQLITE_SNAPSHOT_ROOT}/{PIPELINE_VERSION}"
SQLITE_SNAPSHOT_QUEUE_TOKEN = (
    f"{_OUTPUT_ROOT}/sqlite/sqlite_snapshot_queued.tkn"
    if _OUTPUT_ROOT else
    "sqlite/sqlite_snapshot_queued.tkn"
)

# Phase14 (evidence): evidence/build-gene-report.py writes one self-contained
# annotation report per gene from CONSOLIDATION_MERGED, after scoring and fingerprint.
EVIDENCE_SCRIPTS_DIR = os.path.join(WORKFLOW_DIR, "evidence")
EVIDENCE_PREPARED_DIR = f"{GENOME_PREFIX}evidence/prepared"

# Phase15 (llm): llm/score-genes-llm.py adds an LLM-judged verdict as separate
# columns next to phase11's confidence_score, reading phase14's evidence reports.
LLM_SCRIPTS_DIR = os.path.join(WORKFLOW_DIR, "llm")
LLM_REPORTS_DIR = f"{GENOME_PREFIX}llm/reports"
LLM_SUMMARY = f"{GENOME_PREFIX}llm/llm-summary.tsv"
LLM_COMPUTE_TOKEN = f"{GENOME_PREFIX}llm/llm_compute.tkn"
LLM_TOKEN = f"{GENOME_PREFIX}llm/llm_db.tkn"
# Publication table joining the final table with llm-summary; kept in scoring/.
FINAL_LLM_ANNOTATED_PUBLICATION = f"{GENOME_PREFIX}scoring/FINAL_LLM_labeled-genes-annotated.tsv"

# Phase12 (fingerprint): fingerprint/add-gene-fingerprint.py writes all five
# outputs in one pass and runs after scoring, since full-with-scores needs
# SCORING_CONFIDENCE_FINAL.
FINGERPRINT_SCRIPTS_DIR = os.path.join(WORKFLOW_DIR, "fingerprint")
FINGERPRINT_HASH_PATTERN = f"{GENOME_PREFIX}fingerprint/labeled-genes-fingerprint-hash-pattern.tsv"
FINGERPRINT_HASH_LABEL = f"{GENOME_PREFIX}fingerprint/labeled-genes-fingerprint-hash-label.tsv"
FINGERPRINT_LABEL_PATTERN = f"{GENOME_PREFIX}fingerprint/labeled-genes-fingerprint-label-pattern.tsv"
FINGERPRINT_FULL = f"{GENOME_PREFIX}fingerprint/labeled-genes-fingerprint-full.tsv"
FINGERPRINT_FULL_WITH_SCORES = f"{GENOME_PREFIX}fingerprint/labeled-genes-fingerprint-full-with-scores.tsv"
# Final annotated table: confidence_final plus fingerprint and operon data, in scoring/.
FINAL_ANNOTATED = f"{GENOME_PREFIX}scoring/scored-raw-labeled-genes-final-annotated.tsv"
# Curated publication-ready subset (~43 cols) with full scoring transparency.
FINAL_ANNOTATED_PUBLICATION = f"{GENOME_PREFIX}scoring/FINAL-scored-labeled-genes-annotated.tsv"
# GFF3 annotation file built from FINAL_ANNOTATED_PUBLICATION + rast.gff (no LLM).
ANNOTATION_GFF = f"{GENOME_PREFIX}scoring/annotation.gff3"
FINGERPRINT_COMPUTE_TOKEN = f"{GENOME_PREFIX}fingerprint/fingerprint_compute.tkn"
FINGERPRINT_TOKEN = f"{GENOME_PREFIX}fingerprint/fingerprint_db.tkn"

# Shared cross-genome fingerprint database, one file for all genomes;
# update-fingerprint-database.py serializes updates with fcntl locks.
FINGERPRINT_DATABASE_PATH = _resolve_shared_file(
    'margie_sb.fingerprint_database.path',
    'fingerprint_database.path',
    f'{STORE_ROOT}/fingerprint-database/fingerprint-database.tsv',
    'fingerprint-database.tsv',
)
FINGERPRINT_DATABASE_UPDATED_TOKEN = f"{GENOME_PREFIX}fingerprint/fingerprint_database_updated.tkn"
_FINGERPRINT_DATABASE_DIR = os.path.dirname(FINGERPRINT_DATABASE_PATH)

# Per-operon fingerprint (add-operon-fingerprint.py): groups gene fingerprints
# by operon_id into evidence/label x ordered/composition signals, one row per gene.
OPERON_FINGERPRINT = f"{GENOME_PREFIX}fingerprint/labeled-genes-operon-fingerprint.tsv"
OPERON_FINGERPRINT_DATABASE_EVIDENCE_ORDERED = f"{_FINGERPRINT_DATABASE_DIR}/operon-fingerprint-database-evidence-ordered.tsv"
OPERON_FINGERPRINT_DATABASE_EVIDENCE_COMPOSITION = f"{_FINGERPRINT_DATABASE_DIR}/operon-fingerprint-database-evidence-composition.tsv"
OPERON_FINGERPRINT_DATABASE_LABEL_ORDERED = f"{_FINGERPRINT_DATABASE_DIR}/operon-fingerprint-database-label-ordered.tsv"
OPERON_FINGERPRINT_DATABASE_LABEL_COMPOSITION = f"{_FINGERPRINT_DATABASE_DIR}/operon-fingerprint-database-label-composition.tsv"
OPERON_FINGERPRINT_DATABASE_UPDATED_TOKEN = f"{GENOME_PREFIX}fingerprint/operon_fingerprint_database_updated.tkn"

# Phase13 (synteny/collinearity): containerized ani, aai and closest-organisms.
# ani/aai compare against a persistent cross-run genome pool that each genome's
# .fna and rast.faa are copied into after scoring.
GENOME_POOL_PATH = _resolve_shared_dir(
    'margie_sb.genome_pool.path',
    'genome_pool.path',
    f'{STORE_ROOT}/genome-pool',
)
GENOME_POOL_FNA_DIR = f"{GENOME_POOL_PATH}/fna"
GENOME_POOL_FAA_DIR = f"{GENOME_POOL_PATH}/faa"
GENOME_POOL_TOKEN = f"{GENOME_PREFIX}genome_pool/genome_pool_copy.tkn"

ANI_BATCH_PREFIX = f"{_OUTPUT_ROOT}/original_container_outputs/ani" if _OUTPUT_ROOT else "original_container_outputs/ani"
ANI_BATCH_OUTPUT_DIR = f"{ANI_BATCH_PREFIX}/container_outputs"
ANI_RESULTS = f"{_OUTPUT_ROOT}/ani/ani_results.tsv" if _OUTPUT_ROOT else "ani/ani_results.tsv"
ANI_COMPUTE_TOKEN = f"{_OUTPUT_ROOT}/ani/ani_compute.tkn" if _OUTPUT_ROOT else "ani/ani_compute.tkn"
ANI_TOKEN = f"{_OUTPUT_ROOT}/ani/ani_db.tkn" if _OUTPUT_ROOT else "ani/ani_db.tkn"

AAI_BATCH_PREFIX = f"{_OUTPUT_ROOT}/original_container_outputs/aai" if _OUTPUT_ROOT else "original_container_outputs/aai"
AAI_BATCH_OUTPUT_DIR = f"{AAI_BATCH_PREFIX}/container_outputs"
AAI_RESULTS = f"{_OUTPUT_ROOT}/aai/aai_results.tsv" if _OUTPUT_ROOT else "aai/aai_results.tsv"
AAI_COMPUTE_TOKEN = f"{_OUTPUT_ROOT}/aai/aai_compute.tkn" if _OUTPUT_ROOT else "aai/aai_compute.tkn"
AAI_TOKEN = f"{_OUTPUT_ROOT}/aai/aai_db.tkn" if _OUTPUT_ROOT else "aai/aai_db.tkn"

CLOSEST_BATCH_PREFIX = f"{_OUTPUT_ROOT}/original_container_outputs/closest" if _OUTPUT_ROOT else "original_container_outputs/closest"
CLOSEST_BATCH_OUTPUT_DIR = f"{CLOSEST_BATCH_PREFIX}/container_outputs"
CLOSEST_RESULTS = f"{_OUTPUT_ROOT}/closest/closest_organisms.tsv" if _OUTPUT_ROOT else "closest/closest_organisms.tsv"
CLOSEST_COMPUTE_TOKEN = f"{_OUTPUT_ROOT}/closest/closest_compute.tkn" if _OUTPUT_ROOT else "closest/closest_compute.tkn"
CLOSEST_TOKEN = f"{_OUTPUT_ROOT}/closest/closest_db.tkn" if _OUTPUT_ROOT else "closest/closest_db.tkn"

# Paths for per-genome mauve/synteny against the top-N closest references
# (staged by stage_references.py); no rules in this file use them.
SYNTENY_SCRIPTS_DIR = os.path.join(WORKFLOW_DIR, "synteny")
# Genome name -> raw FASTA path, read by stage_references.py.
GENOME_FASTA_INDEX = f"{_OUTPUT_ROOT}/synteny/genome_fasta_index.tsv" if _OUTPUT_ROOT else "synteny/genome_fasta_index.tsv"
MAUVE_REFERENCES_DIR = f"{GENOME_PREFIX}mauve/references"
MAUVE_RESULTS = f"{GENOME_PREFIX}mauve/conserved_blocks.tsv"
MAUVE_COMPUTE_TOKEN = f"{GENOME_PREFIX}mauve/mauve_compute.tkn"
MAUVE_TOKEN = f"{GENOME_PREFIX}mauve/mauve_db.tkn"

SYNTENY_REFERENCES_DIR = f"{GENOME_PREFIX}synteny/references"
SYNTENY_MERGED_GFF3 = f"{GENOME_PREFIX}synteny/merged_annotation.gff3"
SYNTENY_COMPUTE_TOKEN = f"{GENOME_PREFIX}synteny/synteny_compute.tkn"
SYNTENY_TOKEN = f"{GENOME_PREFIX}synteny/synteny_db.tkn"

# Per genome, phase9 waits for phase4-8 and each later phase depends on the
# previous one's outputs; workflow.py drives phase4_12_one_genome per genome.


def _phase4_8_targets_for_genome(genome):
    """Returns the selected phase4-8 db tokens for one genome, gated by run_<tool>.

    Used by rule all and by phase4_8_one_genome, which workflow.py runs once
    per genome in the order their gene calls finish.
    """
    # An empty genome (no target_genome) returns [] to avoid '//' paths.
    if not genome:
        return []
    targets = []
    for t in PHASE4_TOOLS:
        if rc_bool(f'run_{t}', True, config=config):
            targets.append(PHASE4_TOKENS[t].format(genome=genome))
    if rc_bool('run_interpro', True, config=config):
        targets.extend(
            INTERPRO_PERDB_TOKEN_PATTERN.format(genome=genome, db=db)
            for db in INTERPRO_DB_BASENAMES
        )
    if rc_bool('run_operon', True, config=config):
        targets.append(OPERON_TOKEN.format(genome=genome))
    if rc_bool('run_phobius', True, config=config):
        targets.append(PHOBIUS_TOKEN.format(genome=genome))
    if rc_bool('run_tmbed', True, config=config):
        targets.append(TMBED_TOKEN.format(genome=genome))
    if rc_bool('run_signalp6', True, config=config):
        targets.append(SIGNALP6_TOKEN.format(genome=genome))
    if rc_bool('run_envelope', True, config=config):
        targets.append(ENVELOPE_TOKEN.format(genome=genome))
    if rc_bool('run_deepsig', True, config=config):
        targets.append(DEEPSIG_TOKEN.format(genome=genome))
    if rc_bool('run_psortb', True, config=config):
        targets.append(PSORTB_TOKEN.format(genome=genome))
    if rc_bool('run_signalp4', True, config=config):
        targets.append(SIGNALP4_TOKEN.format(genome=genome))
    return targets


def _phase9_12_targets_for_genome(genome, include_llm=True):
    """Returns the selected phase9-15 targets for one genome, gated by run_<step>.

    Defined before rule all, which calls it at parse time. include_llm=False
    omits the LLM token so the LLM can run separately via rule llm_all.
    """
    if not genome:
        return []
    targets = []
    if rc_bool('run_consolidation', True, config=config):
        targets.append(CONSOLIDATION_TOKEN.format(genome=genome))
    if rc_bool('run_labeling', True, config=config):
        targets.append(LABELING_TOKEN.format(genome=genome))
    if rc_bool('run_scoring', True, config=config):
        targets.append(SCORING_TOKEN.format(genome=genome))
        if rc_bool('run_scoring_archive', True, config=config):
            targets.append(SCORING_ARCHIVE_TOKEN.format(genome=genome))
        if rc_bool('run_final_tables_depot_publish', True, config=config):
            targets.append(FINAL_TABLES_DEPOT_TOKEN.format(genome=genome))
        if rc_bool('run_annotation_gff', False, config=config):
            targets.append(ANNOTATION_GFF.format(genome=genome))
        # Per-organism report figures (default on); never blocks the genome.
        if rc_bool('run_report_figures', True, config=config):
            targets.append(REPORT_FIGURES_ORGANISM_TOKEN.format(genome=genome))
        # Interactive genome viewer (default on); never blocks the genome.
        if rc_bool('run_genome_viewer', True, config=config):
            targets.append(GENOME_VIEWER_TOKEN.format(genome=genome))
        # Full operon atlas (~2-3 min per genome), off by default.
        if rc_bool('run_full_operon_map', _FULL_OPERON_MAP_DEFAULT, config=config):
            targets.append(COMPLETE_OPERON_MAP_TOKEN.format(genome=genome))
    if rc_bool('run_fingerprint', False, config=config):
        targets.append(FINGERPRINT_TOKEN.format(genome=genome))
    if rc_bool('run_fingerprint_database', False, config=config):
        targets.append(FINGERPRINT_DATABASE_UPDATED_TOKEN.format(genome=genome))
    if rc_bool('run_operon_fingerprint', False, config=config):
        targets.append(OPERON_FINGERPRINT_DATABASE_UPDATED_TOKEN.format(genome=genome))
    if rc_bool('run_genome_pool', False, config=config):
        targets.append(GENOME_POOL_TOKEN.format(genome=genome))
    if rc_bool('run_evidence', False, config=config):
        targets.append(EVIDENCE_PREPARED_DIR.format(genome=genome))
    if include_llm and rc_bool('run_llm', False, config=config):
        targets.append(LLM_TOKEN.format(genome=genome))
    return targets


rule all:
    # Requests every selected tool's db token for every genome (run_<tool>, default on).
    input:
        (expand(QUAST_TOKEN, genome=list(GENOMES.keys())) if rc_bool('run_quast', True, config=config) else []),
        (expand(GTDBTK_TOKEN, genome=list(GENOMES.keys())) if rc_bool('run_gtdbtk', True, config=config) else []),
        (expand(RASTTK_TOKEN, genome=list(GENOMES.keys())) if rc_bool('run_rasttk', True, config=config) else []),
        [_phase4_8_targets_for_genome(genome) for genome in GENOMES.keys()],
        [_phase9_12_targets_for_genome(genome) for genome in GENOMES.keys()],
        ([SQLITE_SNAPSHOT_QUEUE_TOKEN]
         if (rc_bool('run_scoring', True, config=config)
             and rc_bool('run_sqlite_snapshot_queue', True, config=config))
         else [])


rule rasttk_all:
    """Stage 1 of workflow.py's sequential run: phase1-3 for every genome.

    run_gtdbtk gates only the GTDB-Tk database load here.
    """
    input:
        expand(RASTTK_TOKEN, genome=list(GENOMES.keys())),
        (expand(GTDBTK_TOKEN, genome=list(GENOMES.keys())) if rc_bool('run_gtdbtk', True, config=config) else [])


rule phase4_8_one_genome:
    """Runs every selected phase4-8 tool for the genome named by --config target_genome."""
    input:
        _phase4_8_targets_for_genome(rc('target_genome', '', config=config))


rule phase4_12_one_genome:
    """Stage 2 of workflow.py's sequential run: phase4-15 for the target_genome."""
    input:
        _phase4_8_targets_for_genome(rc('target_genome', '', config=config)),
        _phase9_12_targets_for_genome(rc('target_genome', '', config=config))


rule phase4_12_one_genome_no_llm:
    """Stage 2 without the LLM token, used when the LLM runs separately in rule llm_all."""
    input:
        _phase4_8_targets_for_genome(rc('target_genome', '', config=config)),
        _phase9_12_targets_for_genome(rc('target_genome', '', config=config), include_llm=False)


rule llm_all:
    """Stage 3: requests the LLM token for every genome; gres=gpu:1 serializes the GPU jobs."""
    input:
        expand(LLM_TOKEN, genome=list(GENOMES.keys())) if rc_bool('run_llm', False, config=config) else []


rule run_consolidation:
    """Phase9: merges every phase4-8 results.tsv into one row per gene, plus a no-stat view.

    Runs detect-columns.py, merge-all-columns.py and filter-no-stat.py on the host.
    """
    input:
        phase4_8=lambda wildcards: _phase4_8_targets_for_genome(wildcards.genome),
        gtdbtk_results=GENOME_INFO
    output:
        detected_columns=CONSOLIDATION_DETECTED_COLUMNS,
        merged=CONSOLIDATION_MERGED,
        manifest=CONSOLIDATION_MANIFEST,
        no_stat=CONSOLIDATION_NO_STAT,
        tkn=CONSOLIDATION_COMPUTE_TOKEN
    threads: rc('consolidation.threads', 1, config=config)
    resources:
        mem_mb=rc('consolidation.mem_mb', 8000, config=config),
        runtime=runtime_min('consolidation.runtime', 30, config=config)
    params:
        genome_dir=lambda wildcards: GENOME_PREFIX.format(genome=wildcards.genome).rstrip('/')
    shell:
        """
        echo "=== MARGIE_SB PHASE 9: CONSOLIDATION ({wildcards.genome}) ==="
        DOMAIN=$(awk -F'\\t' 'NR==1{{for(i=1;i<=NF;i++) if($i=="GTDBTK_domain") c=i}} NR==2{{print $c}}' {input.gtdbtk_results})
        {LOADER_PYTHON} {CONSOLIDATION_SCRIPTS_DIR}/detect-columns.py \
            --input-root {params.genome_dir} \
            --organism-name {wildcards.genome} \
            --output {output.detected_columns}
        {LOADER_PYTHON} {CONSOLIDATION_SCRIPTS_DIR}/merge-all-columns.py \
            --input-root {params.genome_dir} \
            --organism-name {wildcards.genome} \
            --domain "$DOMAIN" \
            --output {output.merged} \
            --manifest {output.manifest}
        {LOADER_PYTHON} {CONSOLIDATION_SCRIPTS_DIR}/filter-no-stat.py \
            --input {output.merged} \
            --output {output.no_stat}
        echo "consolidation complete for {wildcards.genome}" > {output.tkn}
        """


rule run_labeling:
    """Phase10: assigns each gene's canonical_label by trust hierarchy, then adds
    EC consensus, operon info and cluster agreement views (host-side scripts)."""
    input:
        merged=CONSOLIDATION_MERGED,
        consolidation_tkn=CONSOLIDATION_COMPUTE_TOKEN
    output:
        labeled=LABELING_LABELED,
        ec_consensus=LABELING_EC_CONSENSUS,
        operon_info=LABELING_OPERON_INFO,
        cluster_agreement=LABELING_CLUSTER_AGREEMENT,
        tkn=LABELING_COMPUTE_TOKEN
    threads: rc('labeling.threads', 1, config=config)
    resources:
        mem_mb=rc('labeling.mem_mb', 8000, config=config),
        runtime=runtime_min('labeling.runtime', 30, config=config)
    shell:
        """
        echo "=== MARGIE_SB PHASE 10: LABELING ({wildcards.genome}) ==="
        {LOADER_PYTHON} {LABELING_SCRIPTS_DIR}/assign-canonical-label.py \
            --input {input.merged} \
            --output {output.labeled}
        {LOADER_PYTHON} {LABELING_SCRIPTS_DIR}/add-ec-consensus.py \
            --labeled-input {output.labeled} \
            --merged-input {input.merged} \
            --output {output.ec_consensus}
        {LOADER_PYTHON} {LABELING_SCRIPTS_DIR}/add-operon-info.py \
            --labeled-input {output.labeled} \
            --merged-input {input.merged} \
            --output {output.operon_info}
        {LOADER_PYTHON} {LABELING_SCRIPTS_DIR}/add-cluster-agreement.py \
            --labeled-input {output.labeled} \
            --merged-input {input.merged} \
            --output {output.cluster_agreement}
        echo "labeling complete for {wildcards.genome}" > {output.tkn}
        """


rule update_c3_occ_reference_depot:
    """Adds this organism to the shared OCC operon reference before scoring."""
    input:
        labeled=LABELING_LABELED,
        operon_info=LABELING_OPERON_INFO,
        labeling_tkn=LABELING_COMPUTE_TOKEN
    output:
        tkn=SCORING_OCC_REFERENCE_TOKEN
    threads: 1
    resources:
        mem_mb=rc('scoring_occ_update.mem_mb', 4000, config=config),
        runtime=runtime_min('scoring_occ_update.runtime', 30, config=config)
    params:
        reference=C3_REFERENCE_PKL
    shell:
        """
        echo "=== MARGIE_SB PHASE 11: OCC-REFERENCE UPDATE ({wildcards.genome}) ==="
        {LOADER_PYTHON} {SCORING_SCRIPTS_DIR}/update-occ-reference-depot.py \
            --organism {wildcards.genome} \
            --labeled-input {input.labeled} \
            --operon-info-input {input.operon_info} \
            --reference {params.reference} \
            --output-token {output.tkn}
        """


rule run_scoring:
    """Phase11: computes per-gene tiers, C1-C4 and the blended confidence_score.

    C3 is the geometric mean of UniOP pair probability x adjacency reliability
    against the OCC reference. Host-side scripts, run in dependency order.
    """
    input:
        merged=CONSOLIDATION_MERGED,
        labeled=LABELING_LABELED,
        ec_consensus=LABELING_EC_CONSENSUS,
        operon_info=LABELING_OPERON_INFO,
        operon_results=OPERON_RESULTS,
        cluster_agreement=LABELING_CLUSTER_AGREEMENT,
        occ_ref_tkn=SCORING_OCC_REFERENCE_TOKEN,
        labeling_tkn=LABELING_COMPUTE_TOKEN
    output:
        hierarchy_tier=SCORING_HIERARCHY_TIER,
        confidence_tier=SCORING_CONFIDENCE_TIER,
        c1=SCORING_C1,
        c2=SCORING_C2,
        c3=SCORING_C3,
        c4=SCORING_C4,
        confidence_final=SCORING_CONFIDENCE_FINAL,
        final_annotation_with_confidence=SCORING_FINAL_ANNOTATION_WITH_CONFIDENCE,
        tkn=SCORING_COMPUTE_TOKEN
    threads: rc('scoring.threads', 1, config=config)
    resources:
        mem_mb=rc('scoring.mem_mb', 8000, config=config),
        runtime=runtime_min('scoring.runtime', 30, config=config)
    shell:
        """
        echo "=== MARGIE_SB PHASE 11: SCORING ({wildcards.genome}) ==="
        {LOADER_PYTHON} {SCORING_SCRIPTS_DIR}/score-hierarchy-tier.py \
            --labeled-input {input.labeled} \
            --output {output.hierarchy_tier}
        {LOADER_PYTHON} {SCORING_SCRIPTS_DIR}/score-confidence-tier.py \
            --hierarchy-tier-input {output.hierarchy_tier} \
            --ec-consensus-input {input.ec_consensus} \
            --output {output.confidence_tier}
        {LOADER_PYTHON} {SCORING_SCRIPTS_DIR}/score-c1-tool-coverage.py \
            --labeled-input {input.labeled} \
            --cluster-agreement-input {input.cluster_agreement} \
            --output {output.c1}
        {LOADER_PYTHON} {SCORING_SCRIPTS_DIR}/score-c2-operon-probability.py \
            --operon-input {input.operon_info} \
            --operon-results {input.operon_results} \
            --output {output.c2}
        {LOADER_PYTHON} {SCORING_SCRIPTS_DIR}/c3_score_organism.py \
            --operon-info {input.operon_info} \
            --genes-file {input.labeled} \
            --operon-results {input.operon_results} \
            --reference {C3_REFERENCE_PKL} \
            --output {output.c3}
        {LOADER_PYTHON} {SCORING_SCRIPTS_DIR}/score-c4-ec-agreement.py \
            --ec-consensus-input {input.ec_consensus} \
            --confidence-tier-input {output.confidence_tier} \
            --output {output.c4}
        {LOADER_PYTHON} {SCORING_SCRIPTS_DIR}/score-confidence-final.py \
            --c1-input {output.c1} \
            --c2-input {output.c2} \
            --c3-input {output.c3} \
            --c4-input {output.c4} \
            --operon-input {input.operon_info} \
            --merged-input {input.merged} \
            --output {output.confidence_final}
        {LOADER_PYTHON} {SCORING_SCRIPTS_DIR}/make-final-annotation-with-confidence.py \
            --labeled-input {input.labeled} \
            --operon-info-input {input.operon_info} \
            --operon-results-input {input.operon_results} \
            --merged-input {input.merged} \
            --hierarchy-tier-input {output.hierarchy_tier} \
            --confidence-final-input {output.confidence_final} \
            --c4-input {output.c4} \
            --occ-reference {C3_REFERENCE_PKL} \
            --output {output.final_annotation_with_confidence}
        echo "scoring complete for {wildcards.genome}" > {output.tkn}
        """


rule run_fingerprint:
    """Phase12: builds five per-gene fingerprint strings from the labeled table and
    confidence scores (add-gene-fingerprint.py, host-side)."""
    input:
        labeled=LABELING_LABELED,
        confidence_final=SCORING_CONFIDENCE_FINAL,
        scoring_tkn=SCORING_COMPUTE_TOKEN
    output:
        hash_pattern=FINGERPRINT_HASH_PATTERN,
        hash_label=FINGERPRINT_HASH_LABEL,
        label_pattern=FINGERPRINT_LABEL_PATTERN,
        full=FINGERPRINT_FULL,
        full_with_scores=FINGERPRINT_FULL_WITH_SCORES,
        tkn=FINGERPRINT_COMPUTE_TOKEN
    threads: rc('fingerprint.threads', 1, config=config)
    resources:
        mem_mb=rc('fingerprint.mem_mb', 8000, config=config),
        runtime=runtime_min('fingerprint.runtime', 30, config=config)
    shell:
        """
        echo "=== MARGIE_SB PHASE 12: FINGERPRINT ({wildcards.genome}) ==="
        {LOADER_PYTHON} {FINGERPRINT_SCRIPTS_DIR}/add-gene-fingerprint.py \
            --labeled-input {input.labeled} \
            --confidence-final-input {input.confidence_final} \
            --output-hash-pattern {output.hash_pattern} \
            --output-hash-label {output.hash_label} \
            --output-label-pattern {output.label_pattern} \
            --output-full {output.full} \
            --output-full-with-scores {output.full_with_scores}
        echo "fingerprint complete for {wildcards.genome}" > {output.tkn}
        """


rule update_fingerprint_database:
    """Merges this genome's fingerprint hash-label pairs into the shared fingerprint database."""
    input:
        hash_label=FINGERPRINT_HASH_LABEL,
        compute_tkn=FINGERPRINT_COMPUTE_TOKEN
    output:
        tkn=FINGERPRINT_DATABASE_UPDATED_TOKEN
    params:
        db=FINGERPRINT_DATABASE_PATH
    threads: rc('fingerprint_database.threads', 1, config=config)
    resources:
        mem_mb=rc('fingerprint_database.mem_mb', 4000, config=config),
        runtime=runtime_min('fingerprint_database.runtime', 15, config=config)
    shell:
        """
        {LOADER_PYTHON} {FINGERPRINT_SCRIPTS_DIR}/update-fingerprint-database.py \
            --hash-label-input {input.hash_label} \
            --organism {wildcards.genome} \
            --fingerprint-database {params.db}
        echo "fingerprint-database updated for {wildcards.genome}" > {output.tkn}
        """


rule run_operon_fingerprint:
    """Builds operon fingerprints from gene fingerprints grouped by operon_id, then
    the final annotated and publication tables, which need the operon fingerprints."""
    input:
        operon_info=LABELING_OPERON_INFO,
        hash_label=FINGERPRINT_HASH_LABEL,
        fingerprint_full=FINGERPRINT_FULL,
        ec_consensus=LABELING_EC_CONSENSUS,
        confidence_final=SCORING_CONFIDENCE_FINAL,
        fingerprint_tkn=FINGERPRINT_COMPUTE_TOKEN,
        labeling_genes=LABELING_LABELED,
        phobius_top1=PHOBIUS_TOP1
    output:
        operon_fingerprint=OPERON_FINGERPRINT,
        final_annotated=FINAL_ANNOTATED,
        final_annotated_publication=FINAL_ANNOTATED_PUBLICATION,
        tkn=f"{GENOME_PREFIX}fingerprint/operon_fingerprint_compute.tkn"
    params:
        gene_fp_db=FINGERPRINT_DATABASE_PATH,
        operon_fp_db_label_ordered=OPERON_FINGERPRINT_DATABASE_LABEL_ORDERED,
        operon_fp_db_label_composition=OPERON_FINGERPRINT_DATABASE_LABEL_COMPOSITION
    threads: rc('operon_fingerprint.threads', 1, config=config)
    resources:
        mem_mb=rc('operon_fingerprint.mem_mb', 4000, config=config),
        runtime=runtime_min('operon_fingerprint.runtime', 15, config=config)
    shell:
        """
        echo "=== MARGIE_SB PHASE 12: OPERON FINGERPRINT ({wildcards.genome}) ==="
        {LOADER_PYTHON} {FINGERPRINT_SCRIPTS_DIR}/add-operon-fingerprint.py \
            --operon-input {input.operon_info} \
            --hash-label-input {input.hash_label} \
            --output {output.operon_fingerprint}
        {LOADER_PYTHON} {FINGERPRINT_SCRIPTS_DIR}/add-fingerprint-to-final.py \
            --confidence-final-input {input.confidence_final} \
            --fingerprint-hash-label-input {input.hash_label} \
            --operon-fingerprint-input {output.operon_fingerprint} \
            --fingerprint-database {params.gene_fp_db} \
            --operon-fp-label-ordered-database {params.operon_fp_db_label_ordered} \
            --operon-fp-label-composition-database {params.operon_fp_db_label_composition} \
            --ec-consensus-input {input.ec_consensus} \
            --output {output.final_annotated}
        {LOADER_PYTHON} {FINGERPRINT_SCRIPTS_DIR}/make-final-annotated.py \
            --full-evidence-input {output.final_annotated} \
            --labeling-genes-input {input.labeling_genes} \
            --phobius-top1-input {input.phobius_top1} \
            --fingerprint-full-input {input.fingerprint_full} \
            --output {output.final_annotated_publication}
        echo "operon fingerprint complete for {wildcards.genome}" > {output.tkn}
        """


rule update_operon_fingerprint_database:
    """Merges this genome's distinct operons into the four shared operon-fingerprint databases."""
    input:
        operon_fingerprint=OPERON_FINGERPRINT,
        compute_tkn=f"{GENOME_PREFIX}fingerprint/operon_fingerprint_compute.tkn"
    output:
        tkn=OPERON_FINGERPRINT_DATABASE_UPDATED_TOKEN
    params:
        evidence_ordered=OPERON_FINGERPRINT_DATABASE_EVIDENCE_ORDERED,
        evidence_composition=OPERON_FINGERPRINT_DATABASE_EVIDENCE_COMPOSITION,
        label_ordered=OPERON_FINGERPRINT_DATABASE_LABEL_ORDERED,
        label_composition=OPERON_FINGERPRINT_DATABASE_LABEL_COMPOSITION
    threads: rc('operon_fingerprint_database.threads', 1, config=config)
    resources:
        mem_mb=rc('operon_fingerprint_database.mem_mb', 4000, config=config),
        runtime=runtime_min('operon_fingerprint_database.runtime', 15, config=config)
    shell:
        """
        {LOADER_PYTHON} {FINGERPRINT_SCRIPTS_DIR}/update-operon-fingerprint-database.py \
            --operon-fingerprint-input {input.operon_fingerprint} \
            --organism {wildcards.genome} \
            --evidence-ordered-database {params.evidence_ordered} \
            --evidence-composition-database {params.evidence_composition} \
            --label-ordered-database {params.label_ordered} \
            --label-composition-database {params.label_composition}
        echo "operon-fingerprint-database updated for {wildcards.genome}" > {output.tkn}
        """


rule load_consolidation_to_db:
    """Loads the merged and no-stat consolidation tables into main_database (load_to_db.py tsv)."""
    input:
        merged=CONSOLIDATION_MERGED,
        no_stat=CONSOLIDATION_NO_STAT,
        compute_tkn=CONSOLIDATION_COMPUTE_TOKEN,
        fasta=lambda wildcards: GENOMES[wildcards.genome]
    output:
        tkn=CONSOLIDATION_TOKEN
    params:
        db=MAIN_DATABASE,
        script=LOAD_SCRIPT
    shell:
        """
        {LOADER_PYTHON} {params.script} tsv {input.merged} {params.db} consolidation_merged \
            --fasta {input.fasta} --delete-organism {wildcards.genome}
        {LOADER_PYTHON} {params.script} tsv {input.no_stat} {params.db} consolidation_no_stat \
            --fasta {input.fasta} --delete-organism {wildcards.genome}
        echo "consolidation loaded for {wildcards.genome}" > {output.tkn}
        """


rule load_labeling_to_db:
    """Loads each labeling TSV into its own main_database table."""
    input:
        labeled=LABELING_LABELED,
        ec_consensus=LABELING_EC_CONSENSUS,
        operon_info=LABELING_OPERON_INFO,
        cluster_agreement=LABELING_CLUSTER_AGREEMENT,
        compute_tkn=LABELING_COMPUTE_TOKEN,
        fasta=lambda wildcards: GENOMES[wildcards.genome]
    output:
        tkn=LABELING_TOKEN
    params:
        db=MAIN_DATABASE,
        script=LOAD_SCRIPT
    shell:
        """
        {LOADER_PYTHON} {params.script} tsv {input.labeled} {params.db} labeling_labeled \
            --fasta {input.fasta} --delete-organism {wildcards.genome}
        {LOADER_PYTHON} {params.script} tsv {input.ec_consensus} {params.db} labeling_ec_consensus \
            --fasta {input.fasta} --delete-organism {wildcards.genome}
        {LOADER_PYTHON} {params.script} tsv {input.operon_info} {params.db} labeling_operon_info \
            --fasta {input.fasta} --delete-organism {wildcards.genome}
        {LOADER_PYTHON} {params.script} tsv {input.cluster_agreement} {params.db} labeling_cluster_agreement \
            --fasta {input.fasta} --delete-organism {wildcards.genome}
        echo "labeling loaded for {wildcards.genome}" > {output.tkn}
        """


rule load_scoring_to_db:
    """Loads each scoring TSV, including C1-C4, into its own main_database table.

    --delete-organism with --force always replaces the genome's rows, since scores
    change as the OCC reference grows; archive_scoring_to_depot keeps the history.
    """
    input:
        hierarchy_tier=SCORING_HIERARCHY_TIER,
        confidence_tier=SCORING_CONFIDENCE_TIER,
        c1=SCORING_C1,
        c2=SCORING_C2,
        c3=SCORING_C3,
        c4=SCORING_C4,
        confidence_final=SCORING_CONFIDENCE_FINAL,
        compute_tkn=SCORING_COMPUTE_TOKEN,
        fasta=lambda wildcards: GENOMES[wildcards.genome]
    output:
        tkn=SCORING_TOKEN
    params:
        db=MAIN_DATABASE,
        script=LOAD_SCRIPT
    shell:
        """
        {LOADER_PYTHON} {params.script} tsv {input.hierarchy_tier} {params.db} scoring_hierarchy_tier \
            --fasta {input.fasta} --delete-organism {wildcards.genome} --force
        {LOADER_PYTHON} {params.script} tsv {input.confidence_tier} {params.db} scoring_confidence_tier \
            --fasta {input.fasta} --delete-organism {wildcards.genome} --force
        {LOADER_PYTHON} {params.script} tsv {input.c1} {params.db} scoring_c1_tool_coverage \
            --fasta {input.fasta} --delete-organism {wildcards.genome} --force
        {LOADER_PYTHON} {params.script} tsv {input.c2} {params.db} scoring_c2_operon_probability \
            --fasta {input.fasta} --delete-organism {wildcards.genome} --force
        {LOADER_PYTHON} {params.script} tsv {input.c3} {params.db} scoring_c3_operon_context \
            --fasta {input.fasta} --delete-organism {wildcards.genome} --force
        {LOADER_PYTHON} {params.script} tsv {input.c4} {params.db} scoring_c4_ec_agreement \
            --fasta {input.fasta} --delete-organism {wildcards.genome} --force
        {LOADER_PYTHON} {params.script} tsv {input.confidence_final} {params.db} scoring_confidence_final \
            --fasta {input.fasta} --delete-organism {wildcards.genome} --force
        echo "scoring loaded for {wildcards.genome}" > {output.tkn}
        """


rule archive_scoring_to_depot:
    """Copies this genome's confidence-final table into the run's scoring archive folder.

    Gated by run_scoring_archive (default true).
    """
    input:
        confidence_final=SCORING_CONFIDENCE_FINAL,
        compute_tkn=SCORING_COMPUTE_TOKEN
    output:
        tkn=SCORING_ARCHIVE_TOKEN
    params:
        dest_dir=SCORING_ARCHIVE_DIR
    threads: 1
    resources:
        mem_mb=rc('scoring_archive.mem_mb', 1000, config=config),
        runtime=runtime_min('scoring_archive.runtime', 10, config=config)
    shell:
        """
        echo "=== MARGIE_SB: ARCHIVE SCORING TO DEPOT ({wildcards.genome}) ==="
        mkdir -p {params.dest_dir}
        cp {input.confidence_final} {params.dest_dir}/{wildcards.genome}-scored-labeled-genes-confidence-final.tsv
        echo "scoring archived for {wildcards.genome} -> {params.dest_dir}" > {output.tkn}
        """


rule publish_final_annotation_to_depot:
    """Copies the final annotation table to final-tables/<organism>/."""
    input:
        final_with_confidence=SCORING_FINAL_ANNOTATION_WITH_CONFIDENCE,
        scoring_tkn=SCORING_COMPUTE_TOKEN
    output:
        tkn=FINAL_TABLES_DEPOT_TOKEN
    params:
        dest_root=FINAL_TABLES_DEPOT_PATH
    threads: 1
    resources:
        mem_mb=rc('final_tables_depot.mem_mb', 1000, config=config),
        runtime=runtime_min('final_tables_depot.runtime', 10, config=config)
    shell:
        """
        echo "=== MARGIE_SB: PUBLISH FINAL TABLE TO DEPOT ({wildcards.genome}) ==="
        DEST="{params.dest_root}/{wildcards.genome}"
        mkdir -p "$DEST"
        cp {input.final_with_confidence} "$DEST/FINAL_ANNOTATION_WITH_CONFIDENCE.tsv"
        echo "final annotation published for {wildcards.genome} -> $DEST" > {output.tkn}
        """


rule run_report_figures_one_genome:
    """Writes per-organism report figures and TSVs to <genome>/scoring/figures/.

    Errors are logged and ignored so the rule never blocks the genome.
    Gated by run_report_figures (default true).
    """
    input:
        scoring_tkn=SCORING_TOKEN
    output:
        tkn=REPORT_FIGURES_ORGANISM_TOKEN
    params:
        scripts=REPORT_FIGURES_SCRIPTS_DIR,
        outdir=REPORT_FIGURES_ORGANISM_DIR,
        operon_db=REPORT_FIGURES_OPERON_DB,
        run_root=_OUTPUT_ROOT
    threads: 1
    resources:
        mem_mb=rc('report_figures.mem_mb', 8000, config=config),
        runtime=runtime_min('report_figures.runtime', 30, config=config)
    shell:
        """
        echo "=== MARGIE_SB: REPORT FIGURES ({wildcards.genome}) ==="
        {LOADER_PYTHON} {params.scripts}/make_organism_report.py \
            --run-root "{params.run_root}" --organism "{wildcards.genome}" \
            --operon-db "{params.operon_db}" --output-dir "{params.outdir}" \
            || echo "[report_figures] organism figures failed (non-fatal)"
        {LOADER_PYTHON} {params.scripts}/verify_report.py \
            --run-root "{params.run_root}" --organism "{wildcards.genome}" \
            --figures-dir "{params.outdir}" \
            || echo "[report_figures] organism verify reported issues (non-fatal)"
        echo "report figures generated for {wildcards.genome}" > {output.tkn}
        """


rule run_genome_viewer_one_genome:
    """Builds the self-contained HTML genome/operon viewer and the circular map PNG.

    --consolidated is explicit because the generator's default lookup assumes the
    reorganized layout. Errors never fail the rule. Gated by run_genome_viewer.
    """
    input:
        scoring_tkn=SCORING_TOKEN,
        final_with_confidence=SCORING_FINAL_ANNOTATION_WITH_CONFIDENCE,
        merged=CONSOLIDATION_MERGED
    output:
        tkn=GENOME_VIEWER_TOKEN
    params:
        scripts=VIZ_SCRIPTS_DIR,
        html=GENOME_VIEWER_HTML,
        png=GENOME_VIEWER_CIRCULAR_PNG
    threads: 1
    resources:
        mem_mb=rc('genome_viewer.mem_mb', 8000, config=config),
        runtime=runtime_min('genome_viewer.runtime', 30, config=config)
    shell:
        """
        echo "=== MARGIE_SB: GENOME VIEWER ({wildcards.genome}) ==="
        {LOADER_PYTHON} {params.scripts}/gen_genome_viewer.py \
            "{input.final_with_confidence}" "{params.html}" \
            --consolidated "{input.merged}" \
            || echo "[genome_viewer] interactive viewer failed (non-fatal)"
        {LOADER_PYTHON} {params.scripts}/make_circular_genome.py \
            "{input.final_with_confidence}" "{params.png}" \
            || echo "[genome_viewer] circular map failed (non-fatal)"
        echo "genome viewer generated for {wildcards.genome}" > {output.tkn}
        """


rule run_full_operon_map_one_genome:
    """Draws block-arrow maps and gene tables for every multi-gene operon, paginated by size.

    Errors never fail the rule. Gated by run_full_operon_map (default off).
    """
    input:
        scoring_tkn=SCORING_TOKEN
    output:
        tkn=COMPLETE_OPERON_MAP_TOKEN
    params:
        scripts=REPORT_FIGURES_SCRIPTS_DIR,
        outdir=REPORT_FIGURES_ORGANISM_DIR,
        operon_db=REPORT_FIGURES_OPERON_DB,
        run_root=_OUTPUT_ROOT
    threads: 1
    resources:
        mem_mb=rc('full_operon_map.mem_mb', 8000, config=config),
        runtime=runtime_min('full_operon_map.runtime', 90, config=config)
    shell:
        """
        echo "=== MARGIE_SB: FULL OPERON MAP ({wildcards.genome}) ==="
        {LOADER_PYTHON} {params.scripts}/make_complete_operon_diagrams.py \
            --run-root "{params.run_root}" --organism "{wildcards.genome}" \
            --operon-db "{params.operon_db}" --output-dir "{params.outdir}" \
            || echo "[full_operon_map] atlas failed (non-fatal)"
        echo "full operon map generated for {wildcards.genome}" > {output.tkn}
        """


rule run_report_figures_global:
    """Writes pangenome report figures to <run>/scoring/figures/global/ after all genomes score.

    Errors never fail the rule; workflow.py runs it as a separate finalize step.
    """
    input:
        scoring_tokens=expand(SCORING_TOKEN, genome=list(GENOMES.keys()))
    output:
        tkn=REPORT_FIGURES_GLOBAL_TOKEN
    params:
        scripts=REPORT_FIGURES_SCRIPTS_DIR,
        outdir=REPORT_FIGURES_GLOBAL_DIR,
        operon_db=REPORT_FIGURES_OPERON_DB,
        run_root=_OUTPUT_ROOT
    threads: 1
    resources:
        mem_mb=rc('report_figures.global_mem_mb', 16000, config=config),
        runtime=runtime_min('report_figures.global_runtime', 45, config=config)
    shell:
        """
        echo "=== MARGIE_SB: REPORT FIGURES (global / pangenome) ==="
        {LOADER_PYTHON} {params.scripts}/make_global_report.py \
            --run-root "{params.run_root}" --operon-db "{params.operon_db}" \
            --output-dir "{params.outdir}" \
            || echo "[report_figures] global figures failed (non-fatal)"
        {LOADER_PYTHON} {params.scripts}/verify_report.py \
            --run-root "{params.run_root}" --figures-dir "{params.outdir}" \
            || echo "[report_figures] global verify reported issues (non-fatal)"
        echo "global report figures generated" > {output.tkn}
        """


rule queue_sqlite_backup_snapshot:
    """Submits a background SLURM job (sbatch) that copies margie.db into a pipeline-version snapshot."""
    input:
        scoring_tokens=expand(SCORING_TOKEN, genome=list(GENOMES.keys())),
        depot_tokens=(expand(FINAL_TABLES_DEPOT_TOKEN, genome=list(GENOMES.keys()))
                      if rc_bool('run_final_tables_depot_publish', True, config=config) else [])
    output:
        tkn=SQLITE_SNAPSHOT_QUEUE_TOKEN
    params:
        source_db=MAIN_DATABASE,
        dest_dir=SQLITE_SNAPSHOT_VERSION_DIR,
        account=rc('sqlite_backup.account', '', config=config),
        partition=rc('sqlite_backup.partition', '', config=config),
        time_limit=rc('sqlite_backup.time', '01:00:00', config=config)
    threads: 1
    resources:
        mem_mb=rc('sqlite_backup.mem_mb', 1000, config=config),
        runtime=runtime_min('sqlite_backup.runtime', 10, config=config)
    shell:
        """
        echo "=== MARGIE_SB: QUEUE SQLITE BACKUP SNAPSHOT ==="
        mkdir -p $(dirname {output.tkn})
        mkdir -p {params.dest_dir}
        SBATCH_ACCOUNT=""
        SBATCH_PARTITION=""
        if [ -n "{params.account}" ]; then
            SBATCH_ACCOUNT="--account={params.account}"
        fi
        if [ -n "{params.partition}" ]; then
            SBATCH_PARTITION="--partition={params.partition}"
        fi
        DEST_DB="{params.dest_dir}/margie.db"
        JOB_ID=$(sbatch --parsable $SBATCH_ACCOUNT $SBATCH_PARTITION --time={params.time_limit} \
            --job-name=margie-sqlite-snapshot \
            --output={params.dest_dir}/sqlite-backup-%j.out \
            --error={params.dest_dir}/sqlite-backup-%j.err \
            --wrap="set -euo pipefail; mkdir -p {params.dest_dir}; cp {params.source_db} $DEST_DB")
        {{
            echo "will copy the sqlite database to backup location"
            echo "sqlite backup job submitted: $JOB_ID"
            echo "source_db={params.source_db}"
            echo "dest_db=$DEST_DB"
            echo "workflow is complete now."
        }} > {output.tkn}
        cat {output.tkn}
        """






rule make_annotation_gff:
    """Builds scoring/annotation.gff3 (one CDS per gene with labels, fingerprint and
    confidence as attributes) from the publication table and rast.gff contigs."""
    input:
        final=FINAL_ANNOTATED_PUBLICATION,
        rast_gff=RASTTK_GFF,
        scoring_tkn=SCORING_TOKEN
    output:
        gff=ANNOTATION_GFF
    shell:
        """
        echo "=== MARGIE_SB: ANNOTATION GFF ({wildcards.genome}) ==="
        {LOADER_PYTHON} {SCORING_SCRIPTS_DIR}/make-gff.py \
            --final {input.final} \
            --rast-gff {input.rast_gff} \
            --output {output.gff}
        """


rule build_gene_evidence_report:
    """Phase14: writes one annotation report per gene with all database evidence,
    C1-C4 and gene/operon fingerprints with cross-genome frequencies (host-side).

    Kept separate from run_llm so the reports stand on their own.
    """
    input:
        merged=CONSOLIDATION_MERGED,
        confidence_final=SCORING_CONFIDENCE_FINAL,
        fingerprint_full=FINGERPRINT_FULL,
        operon_fingerprint=OPERON_FINGERPRINT,
        fingerprint_tkn=FINGERPRINT_COMPUTE_TOKEN
    output:
        prepared_dir=directory(EVIDENCE_PREPARED_DIR)
    threads: rc('evidence.threads', 1, config=config)
    resources:
        mem_mb=rc('evidence.mem_mb', 8000, config=config),
        runtime=runtime_min('evidence.runtime', 60, config=config)
    params:
        # Shared databases are params, not inputs, so updates by other genomes
        # do not trigger reruns.
        gene_fp_db=FINGERPRINT_DATABASE_PATH,
        operon_fp_db_evidence_ordered=OPERON_FINGERPRINT_DATABASE_EVIDENCE_ORDERED,
        operon_fp_db_evidence_composition=OPERON_FINGERPRINT_DATABASE_EVIDENCE_COMPOSITION,
        operon_fp_db_label_ordered=OPERON_FINGERPRINT_DATABASE_LABEL_ORDERED,
        operon_fp_db_label_composition=OPERON_FINGERPRINT_DATABASE_LABEL_COMPOSITION
    shell:
        """
        echo "=== MARGIE_SB PHASE 14: EVIDENCE ({wildcards.genome}) ==="
        {LOADER_PYTHON} {EVIDENCE_SCRIPTS_DIR}/build-gene-report.py \
            --consolidated {input.merged} \
            --confidence-final {input.confidence_final} \
            --fingerprint-full {input.fingerprint_full} \
            --operon-fingerprint {input.operon_fingerprint} \
            --fingerprint-database {params.gene_fp_db} \
            --operon-fingerprint-database-evidence-ordered {params.operon_fp_db_evidence_ordered} \
            --operon-fingerprint-database-evidence-composition {params.operon_fp_db_evidence_composition} \
            --operon-fingerprint-database-label-ordered {params.operon_fp_db_label_ordered} \
            --operon-fingerprint-database-label-composition {params.operon_fp_db_label_composition} \
            --organism-name {wildcards.genome} \
            --output-dir {output.prepared_dir}
        """


rule run_llm:
    """Phase15 (GPU, llm.sif): asks the model for one verdict per gene on whether its
    confidence_score and canonical_label fit the phase14 evidence report."""
    input:
        prepared_dir=EVIDENCE_PREPARED_DIR
    output:
        reports_dir=directory(LLM_REPORTS_DIR),
        summary=LLM_SUMMARY,
        tkn=LLM_COMPUTE_TOKEN
    threads: rc('llm.threads', 10, config=config)
    resources:
        mem_mb=rc('llm.mem_mb', 32000, config=config),
        runtime=runtime_min('llm.runtime', 240, config=config),
        slurm_partition=rc('llm.partition', 'gpu', config=config),
        gres=rc('llm.gres', 'gpu:1', config=config)
    params:
        # Model weights path from llm.model_path (not db_path, whose layout differs).
        model=rc('llm.model_path', '', config=config),
        batch_size=rc('llm.batch_size', 5, config=config),
        max_tokens=rc('llm.max_tokens', 800, config=config)
    container: sif_path('llm.sif', config=config, workflow_id='margie_sb')
    shell:
        """
        echo "=== MARGIE_SB PHASE 15: LLM (review, {wildcards.genome}) ==="
        export PYTORCH_HIP_ALLOC_CONF=expandable_segments:True
        python3 {LLM_SCRIPTS_DIR}/score-genes-llm.py \
            --trained-model {params.model} \
            --prepared-dir {input.prepared_dir} \
            --reports-dir {output.reports_dir} \
            --summary {output.summary} \
            --batch-size {params.batch_size} \
            --max-tokens {params.max_tokens}
        echo "llm analysis complete for {wildcards.genome}" > {output.tkn}
        """


rule load_llm_to_db:
    """Loads the LLM summary into main_database and builds the LLM publication table.

    Per-gene text reports stay on disk and are not loaded.
    """
    input:
        summary=LLM_SUMMARY,
        compute_tkn=LLM_COMPUTE_TOKEN,
        full_annotated=FINAL_ANNOTATED
    output:
        final_llm_publication=FINAL_LLM_ANNOTATED_PUBLICATION,
        tkn=LLM_TOKEN
    params:
        db=MAIN_DATABASE,
        script=LOAD_SCRIPT
    shell:
        """
        {LOADER_PYTHON} {params.script} tsv {input.summary} {params.db} llm_summary
        {LOADER_PYTHON} {LLM_SCRIPTS_DIR}/make-final-llm-annotated.py \
            --full-annotated-input {input.full_annotated} \
            --llm-summary-input {input.summary} \
            --output {output.final_llm_publication}
        echo "llm loaded for {wildcards.genome}" > {output.tkn}
        """


rule load_fingerprint_to_db:
    """Loads each fingerprint TSV into its own main_database table.

    --delete-organism replaces a genome's rows; full_with_scores also passes
    --force because it embeds the changing confidence score.
    """
    input:
        hash_pattern=FINGERPRINT_HASH_PATTERN,
        hash_label=FINGERPRINT_HASH_LABEL,
        label_pattern=FINGERPRINT_LABEL_PATTERN,
        full=FINGERPRINT_FULL,
        full_with_scores=FINGERPRINT_FULL_WITH_SCORES,
        compute_tkn=FINGERPRINT_COMPUTE_TOKEN
    output:
        tkn=FINGERPRINT_TOKEN
    params:
        db=MAIN_DATABASE,
        script=LOAD_SCRIPT
    shell:
        """
        {LOADER_PYTHON} {params.script} tsv {input.hash_pattern} {params.db} fingerprint_hash_pattern --delete-organism {wildcards.genome}
        {LOADER_PYTHON} {params.script} tsv {input.hash_label} {params.db} fingerprint_hash_label --delete-organism {wildcards.genome}
        {LOADER_PYTHON} {params.script} tsv {input.label_pattern} {params.db} fingerprint_label_pattern --delete-organism {wildcards.genome}
        {LOADER_PYTHON} {params.script} tsv {input.full} {params.db} fingerprint_full --delete-organism {wildcards.genome}
        {LOADER_PYTHON} {params.script} tsv {input.full_with_scores} {params.db} fingerprint_full_with_scores --delete-organism {wildcards.genome} --force
        echo "fingerprint loaded for {wildcards.genome}" > {output.tkn}
        """


rule copy_to_genome_pool:
    """Copies this genome's .fna and rast.faa into the shared genome pool after scoring.

    Each genome writes its own {genome}.fna/.faa, so no locking is needed.
    """
    input:
        fna=lambda wildcards: GENOMES[wildcards.genome],
        faa=RASTTK_FAA,
        scoring_tkn=SCORING_COMPUTE_TOKEN
    output:
        tkn=GENOME_POOL_TOKEN
    threads: 1
    resources:
        mem_mb=rc('genome_pool.mem_mb', 1000, config=config),
        runtime=runtime_min('genome_pool.runtime', 10, config=config)
    params:
        fna_dir=GENOME_POOL_FNA_DIR,
        faa_dir=GENOME_POOL_FAA_DIR
    shell:
        """
        mkdir -p {params.fna_dir} {params.faa_dir}
        cp -f {input.fna} {params.fna_dir}/{wildcards.genome}.fna
        cp -f {input.faa} {params.faa_dir}/{wildcards.genome}.faa
        echo "copied {wildcards.genome} to genome pool" > {output.tkn}
        """


rule run_ani_batch:
    """Phase13: skani all-vs-all nucleotide identity over the whole genome pool,
    after this run's genomes are copied in."""
    input:
        pool_tkns=expand(GENOME_POOL_TOKEN, genome=list(GENOMES.keys()))
    output:
        results=ANI_RESULTS,
        tkn=ANI_COMPUTE_TOKEN
    group: "ani"
    threads: rc('ani.threads', 8, config=config)
    resources:
        mem_mb=rc('ani.mem_mb', 8000, config=config),
        runtime=runtime_min('ani.runtime', 60, config=config)
    params:
        pool_dir=GENOME_POOL_FNA_DIR,
        output_dir=rc('ani.output_dir', ANI_BATCH_OUTPUT_DIR, config=config),
        genome_count=len(GENOMES),
        min_af=rc('ani.min_af', '0.15', config=config),
        kmer=rc('ani.kmer', '15', config=config)
    container: sif_path('ani.sif', config=config, workflow_id='margie_sb')
    shell:
        """
        echo "=== MARGIE_SB PHASE 13: ANI (pool, this run contributed {params.genome_count} genomes) ==="
        /usr/local/bin/run -i {params.pool_dir} -o {params.output_dir} -t {threads} \
            --min-af {params.min_af} --kmer {params.kmer} --organism-name margie_genome_pool
        cp {params.output_dir}/ani/processed/ani_results.tsv {output.results}
        echo "ani complete" > {output.tkn}
        """


rule load_ani_to_db:
    """Loads ANI results into main_database."""
    input:
        results=ANI_RESULTS,
        compute_tkn=ANI_COMPUTE_TOKEN
    output:
        tkn=ANI_TOKEN
    group: "ani"
    params:
        db=MAIN_DATABASE,
        script=LOAD_SCRIPT
    shell:
        """
        {LOADER_PYTHON} {params.script} tsv {input.results} {params.db} ani
        echo "ani loaded" > {output.tkn}
        """


rule run_aai_batch:
    """Phase13: CompareM/DIAMOND all-vs-all amino acid identity over the genome pool's .faa files."""
    input:
        pool_tkns=expand(GENOME_POOL_TOKEN, genome=list(GENOMES.keys()))
    output:
        results=AAI_RESULTS,
        tkn=AAI_COMPUTE_TOKEN
    group: "aai"
    threads: rc('aai.threads', 8, config=config)
    resources:
        mem_mb=rc('aai.mem_mb', 8000, config=config),
        runtime=runtime_min('aai.runtime', 60, config=config)
    params:
        pool_dir=GENOME_POOL_FAA_DIR,
        output_dir=rc('aai.output_dir', AAI_BATCH_OUTPUT_DIR, config=config),
        genome_count=len(GENOMES),
        min_identity=rc('aai.min_identity', '30.0', config=config),
        min_aln_len=rc('aai.min_aln_len', '70.0', config=config),
        evalue=rc('aai.evalue', '0.001', config=config)
    container: sif_path('aai.sif', config=config, workflow_id='margie_sb')
    shell:
        """
        echo "=== MARGIE_SB PHASE 13: AAI (pool, this run contributed {params.genome_count} genomes) ==="
        /usr/local/bin/run -i {params.pool_dir} -o {params.output_dir} -t {threads} \
            --min-identity {params.min_identity} --min-aln-len {params.min_aln_len} \
            --evalue {params.evalue} --organism-name margie_genome_pool
        cp {params.output_dir}/aai/processed/aai_results.tsv {output.results}
        echo "aai complete" > {output.tkn}
        """


rule load_aai_to_db:
    """Loads AAI results into main_database."""
    input:
        results=AAI_RESULTS,
        compute_tkn=AAI_COMPUTE_TOKEN
    output:
        tkn=AAI_TOKEN
    group: "aai"
    params:
        db=MAIN_DATABASE,
        script=LOAD_SCRIPT
    shell:
        """
        {LOADER_PYTHON} {params.script} tsv {input.results} {params.db} aai
        echo "aai loaded" > {output.tkn}
        """


rule run_closest_organisms_batch:
    """Phase13: ranks each genome's top-N closest relatives from ANI, falling back
    to AAI, in one batch call."""
    input:
        ani=ANI_RESULTS,
        aai=AAI_RESULTS,
        ani_tkn=ANI_COMPUTE_TOKEN,
        aai_tkn=AAI_COMPUTE_TOKEN
    output:
        results=CLOSEST_RESULTS,
        tkn=CLOSEST_COMPUTE_TOKEN
    group: "closest"
    threads: rc('closest.threads', 2, config=config)
    resources:
        mem_mb=rc('closest.mem_mb', 4000, config=config),
        runtime=runtime_min('closest.runtime', 30, config=config)
    params:
        output_dir=rc('closest.output_dir', CLOSEST_BATCH_OUTPUT_DIR, config=config),
        top_n=rc('closest.top_n', '5', config=config),
        min_identity=rc('closest.min_identity', '30.0', config=config),
        gap_scaling=rc('closest.gap_scaling', '100', config=config)
    container: sif_path('closest.sif', config=config, workflow_id='margie_sb')
    shell:
        """
        echo "=== MARGIE_SB PHASE 13: CLOSEST-ORGANISMS ==="
        /usr/local/bin/run -o {params.output_dir} --ani {input.ani} --aai {input.aai} \
            --top-n {params.top_n} --min-identity {params.min_identity} \
            --gap-scaling {params.gap_scaling} --collection-name margie_collection
        cp {params.output_dir}/closest/processed/closest_organisms.tsv {output.results}
        echo "closest-organisms complete" > {output.tkn}
        """


rule load_closest_to_db:
    """Loads closest-organisms results into main_database."""
    input:
        results=CLOSEST_RESULTS,
        compute_tkn=CLOSEST_COMPUTE_TOKEN
    output:
        tkn=CLOSEST_TOKEN
    group: "closest"
    params:
        db=MAIN_DATABASE,
        script=LOAD_SCRIPT
    shell:
        """
        {LOADER_PYTHON} {params.script} tsv {input.results} {params.db} closest_organisms
        echo "closest-organisms loaded" > {output.tkn}
        """


rule run_quast_batch:
    """Phase1: runs QUAST once over a staged directory of all genomes."""
    input:
        list(GENOMES.values())
    output:
        done=QUAST_BATCH_DONE
    group: "quast"
    threads: rc('quast.threads', 8, config=config)
    resources:
        mem_mb=rc('quast.mem_mb', 4000, config=config),
        runtime=runtime_min('quast.runtime', 120, config=config)
    params:
        stage_dir=QUAST_BATCH_STAGE_DIR,
        output_dir=rc('quast.output_dir', QUAST_BATCH_OUTPUT_DIR, config=config),
        genome_count=len(GENOMES)
    container: sif_path('quast.sif', config=config, workflow_id='margie_sb')
    shell:
        """
        echo "=== MARGIE_SB PHASE 1: QUAST (batch {params.genome_count} genomes) ==="
        rm -rf {params.stage_dir}
        mkdir -p {params.stage_dir}

        # Stage all discovered genomes into one folder for a single QUAST run.
        i=0
        for src in {input}; do
            base=$(basename "$src")
            dest="{params.stage_dir}/$base"
            if [[ -e "$dest" ]]; then
                stem="${{base%.*}}"
                ext="${{base##*.}}"
                dest="{params.stage_dir}/${{stem}}_dup${{i}}.${{ext}}"
            fi
            cp -f "$src" "$dest"
            i=$((i + 1))
        done

        /usr/local/bin/run -i {params.stage_dir} -o {params.output_dir} -t {threads}
        touch {output.done}
        """


rule split_quast_batch_per_genome:
    """Copies each genome's batched QUAST output to its per-genome path."""
    input:
        done=QUAST_BATCH_DONE
    output:
        results=expand(QUAST_RESULTS, genome=list(GENOMES.keys()))
    params:
        output_dir=rc('quast.output_dir', QUAST_BATCH_OUTPUT_DIR, config=config)
    run:
        from pathlib import Path

        for out_file in output.results:
            target = Path(out_file)
            genome = target.parts[-3]  # .../<genome>/quast/quast.tsv
            source = Path(params.output_dir) / genome / 'processed' / 'quast.tsv'
            if not source.exists():
                raise ValueError(f"Missing QUAST processed output for genome '{genome}': {source}")
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text(source.read_text())


rule load_quast_to_db:
    """Loads QUAST results into main_database."""
    input:
        results=QUAST_RESULTS
    output:
        tkn=QUAST_TOKEN
    group: "quast"
    params:
        db=MAIN_DATABASE,
        script=LOAD_SCRIPT
    shell:
        """
        {LOADER_PYTHON} {params.script} tsv {input.results} {params.db} quast --token {output.tkn}
        """


rule run_gtdbtk_batch:
    """Phase2: runs GTDB-Tk classify_wf once over all genomes and writes the combined
    results and translation-table files."""
    input:
        list(GENOMES.values())
    output:
        results=GTDBTK_BATCH_RESULTS,
        translation_table=GTDBTK_BATCH_TRANSLATION_TABLE,
        done=GTDBTK_BATCH_DONE
    threads: rc('margie_sb.gtdbtk.threads', rc('gtdbtk.threads', 64, config=config), config=config)
    resources:
        mem_mb=rc('margie_sb.gtdbtk.mem_mb',
                  rc('gtdbtk.mem_mb', 460000, config=config),
                  config=config),
        runtime=runtime_min('margie_sb.gtdbtk.runtime',
                   rc('gtdbtk.runtime', 240, config=config),
                   config=config),
        # Tool key, then phase2-wide partition, then the legacy top-level key.
        slurm_partition=rc('margie_sb.gtdbtk.partition',
                   rc('margie_sb.phase2.partition',
                      rc('gtdbtk.partition', 'highmem', config=config),
                      config=config),
                           config=config)
    params:
        stage_dir=GTDBTK_BATCH_STAGE_DIR,
        output_dir=rc('gtdbtk.output_dir', GTDBTK_BATCH_OUTPUT_DIR, config=config),
        db=db_path('gtdbtk', config=config, workflow_id='margie_sb'),
        genome_count=len(GENOMES)
    container: sif_path('gtdbtk.sif', config=config, workflow_id='margie_sb')
    shell:
        """
        echo "=== MARGIE_SB PHASE 2: GTDBTK (batch {params.genome_count} genomes) ==="
        rm -rf {params.stage_dir}
        mkdir -p {params.stage_dir}

        # Stage all discovered genomes into one folder for a single GTDB-Tk
        # classify_wf call. Preserve duplicate basenames by suffixing.
        #
        # Every copy is renamed to .fna. The container passes --extension fna
        # to classify_wf, so anything staged under another suffix is accepted
        # by the container's own staging and then silently ignored by GTDB-Tk
        # -- discover_genomes takes .fasta/.fa/.fna (GENOME_EXTENSIONS), so a
        # single .fasta in the pool left the batch one genome short and only
        # surfaced much later, as a missing row in
        # split_gtdbtk_batch_per_genome. The stem (basename minus the final
        # suffix, and minus .gz first) is exactly the genome key downstream
        # rules use, and is also the name GTDB-Tk reports, so renaming the
        # extension keeps every contract intact. Gzipped inputs are expanded
        # rather than copied, since a .fna holding gzip bytes would parse as
        # an empty genome.
        i=0
        for src in {input}; do
            base=$(basename "$src")
            stem="${{base%.gz}}"
            stem="${{stem%.*}}"
            dest="{params.stage_dir}/${{stem}}.fna"
            if [[ -e "$dest" ]]; then
                dest="{params.stage_dir}/${{stem}}_dup${{i}}.fna"
            fi
            case "$base" in
                *.gz) zcat "$src" > "$dest" ;;
                *)    cp -f "$src" "$dest" ;;
            esac
            i=$((i + 1))
        done

        # Always force full species placement (skips GTDB-Tk's ANI-only fast
        # path) so identify/MSA runs for every genome -- ANI-only genomes
        # never get an identify-stage translation table, which silently
        # produced wrong genetic codes downstream for genomes needing a
        # non-standard table (e.g. Mycoplasmatales' code 4).
        #
        # GTDBTK_PPLACER_CPUS: the container defaults pplacer to 1 CPU
        # regardless of --cpus/-t, since its memory use on GTDB-Tk's
        # bacterial reference tree scales poorly with thread count. Set to 16
        # as a cautious bump; check peak memory via
        # `seff <slurm_job_id>` (MaxRSS) before raising further.
        GTDBTK_PLACE_SPECIES=1 GTDBTK_PPLACER_CPUS=16 /usr/local/bin/run -i {params.stage_dir} -o {params.output_dir} -d {params.db} -t {threads} --collection-name margie_sb_batch

        cp $(find {params.output_dir} -name "gtdbtk_results.tsv") {output.results}
        cp $(find {params.output_dir} -name "gtdbtk.translation_table_summary.tsv") {output.translation_table}
        touch {output.done}
        """


rule split_gtdbtk_batch_per_genome:
    """Writes one genome's gtdbtk_results.tsv and translation_table.tsv from the batch files.

    Per-genome, so a genome restored from cache skips the split and starts RASTtk at once.
    """
    input:
        results=GTDBTK_BATCH_RESULTS,
        translation_table=GTDBTK_BATCH_TRANSLATION_TABLE,
        done=GTDBTK_BATCH_DONE
    output:
        results=GTDBTK_RESULTS,
        translation_table=GTDBTK_TRANSLATION_TABLE,
        tkn=GTDBTK_COMPUTE_TOKEN
    run:
        import csv
        from pathlib import Path

        genome = wildcards.genome

        with open(input.results, newline='') as fh:
            reader = csv.DictReader(fh, delimiter='\t')
            headers = reader.fieldnames or []
            if 'genome' not in headers:
                raise ValueError(f"Expected 'genome' column in {input.results}, got {headers}")
            rows_by_genome = {}
            for row in reader:
                g = row.get('genome', '')
                if g:
                    rows_by_genome[g] = row

        row = rows_by_genome.get(genome)
        if row is None:
            raise ValueError(f"Missing GTDB-Tk row for genome '{genome}' in {input.results}")
        target = Path(output.results)
        target.parent.mkdir(parents=True, exist_ok=True)
        with target.open('w', newline='') as out_fh:
            # lineterminator='\n': csv defaults to CRLF, which breaks awk header matches downstream
            writer = csv.DictWriter(out_fh, fieldnames=list(row.keys()), delimiter='\t', lineterminator='\n')
            writer.writeheader()
            writer.writerow(row)

        rows_by_genome = {}
        with open(input.translation_table, newline='') as fh:
            first = fh.readline().strip()
            rest = fh.read()

        # Typical file has headers (user_genome/genome). Some GTDB-Tk runs
        # emit a compact two-column file without headers: <genome>\t<table>.
        header_like = ('user_genome' in first) or ('genome' in first)
        if header_like:
            import io
            buf = io.StringIO(first + "\n" + rest)
            reader = csv.DictReader(buf, delimiter='\t')
            headers = reader.fieldnames or []
            genome_col = 'user_genome' if 'user_genome' in headers else ('genome' if 'genome' in headers else headers[0] if headers else None)
            if not genome_col:
                raise ValueError(f"Unable to detect genome column in {input.translation_table}")
            for row in reader:
                g = row.get(genome_col, '')
                if g:
                    rows_by_genome[g] = row
        else:
            lines = [first] + ([ln for ln in rest.splitlines() if ln.strip()] if rest else [])
            for ln in lines:
                parts = ln.split('\t')
                if len(parts) >= 2:
                    rows_by_genome[parts[0]] = {
                        'user_genome': parts[0],
                        'translation_table': parts[1],
                    }

        row = rows_by_genome.get(genome)
        if row is None:
            # Full species placement gives every genome a row, so a gap is an error.
            raise ValueError(f"Missing GTDB-Tk translation table row for genome '{genome}' in {input.translation_table}")
        target = Path(output.translation_table)
        target.parent.mkdir(parents=True, exist_ok=True)
        with target.open('w', newline='') as out_fh:
            # lineterminator='\n': csv defaults to CRLF, which breaks awk header matches downstream
            writer = csv.DictWriter(out_fh, fieldnames=list(row.keys()), delimiter='\t', lineterminator='\n')
            writer.writeheader()
            writer.writerow(row)

        Path(output.tkn).write_text(f"gtdbtk split complete for {genome}\n")


rule load_gtdbtk_to_db:
    """Loads GTDB-Tk results into main_database.

    Not grouped with the GTDB-Tk run, which needs the highmem partition.
    """
    input:
        results=GTDBTK_RESULTS
    output:
        tkn=GTDBTK_TOKEN
    params:
        db=MAIN_DATABASE,
        script=LOAD_SCRIPT
    shell:
        """
        {LOADER_PYTHON} {params.script} tsv {input.results} {params.db} gtdbtk --token {output.tkn}
        """


# Runs locally: it only writes one small file.
localrules: resolve_genome_info


rule resolve_genome_info:
    """Writes GENOME_INFO: the genome's domain, genetic code and gene caller.

    Values from margie_sb.genome_info win over GTDB-Tk's; with GTDB-Tk off it has
    no GTDB-Tk inputs. An unknown domain is written as Unknown, an unknown code empty.
    """
    input:
        unpack(lambda wildcards: {'gtdbtk_results': GTDBTK_RESULTS.format(genome=wildcards.genome),
                                  'translation_table': GTDBTK_TRANSLATION_TABLE.format(genome=wildcards.genome)}
               if RUN_GTDBTK else {})
    output:
        info=GENOME_INFO
    run:
        import csv
        from pathlib import Path

        call = GENOME_CALLS[wildcards.genome]
        domain, code, source = call['domain'], call['genetic_code'], call['source']

        def first_value(path, column):
            with open(path, newline='') as fh:
                for row in csv.DictReader(fh, delimiter='\t'):
                    value = (row.get(column) or '').strip()
                    if value:
                        return value
            return ''

        if RUN_GTDBTK:
            # The config only overrides what it actually names.
            domain = domain or first_value(input.gtdbtk_results, 'GTDBTK_domain')
            code = code or first_value(input.translation_table, 'translation_table')
            if not call['domain'] and not call['genetic_code']:
                source = 'gtdbtk'

        target = Path(output.info)
        target.parent.mkdir(parents=True, exist_ok=True)
        with target.open('w', newline='') as fh:
            writer = csv.writer(fh, delimiter='\t', lineterminator='\n')
            writer.writerow(['genome', 'GTDBTK_domain', 'translation_table', 'gene_caller', 'source'])
            writer.writerow([wildcards.genome, domain or 'Unknown', code, call['gene_caller'], source])


rule run_rasttk:
    """Phase3: RASTtk/BV-BRC gene calling and annotation using the genome's domain
    and genetic code from GENOME_INFO.

    BV-BRC calls are serialized by an atomic mkdir lock (RASTTK_BVBRC_LOCK);
    a lock older than this rule's runtime is treated as stale and taken over.
    """
    input:
        fasta=lambda wildcards: GENOMES[wildcards.genome],
        gtdbtk_results=GENOME_INFO,
        translation_table=GENOME_INFO
    output:
        results=RASTTK_RESULTS,
        faa=RASTTK_FAA,
        gff=RASTTK_GFF,
        tkn=RASTTK_COMPUTE_TOKEN
    # RASTtk genomes only; run_prodigal makes the same files for the rest.
    wildcard_constraints:
        genome=_one_of(RASTTK_GENOMES)
    threads: rc('rasttk.threads', 8, config=config)
    resources:
        mem_mb=rc('rasttk.mem_mb', 8000, config=config),
        runtime=runtime_min('rasttk.runtime', 120, config=config)
    params:
        output_dir=lambda wildcards: rc('rasttk.output_dir', f"{CONTAINER_OUTPUTS_PREFIX}rasttk".format(genome=wildcards.genome), config=config),
        db=db_path('rasttk', config=config, workflow_id='margie_sb'),
        lock=RASTTK_BVBRC_LOCK
    container: sif_path('rasttk.sif', config=config, workflow_id='margie_sb')
    shell:
        """
        echo "=== MARGIE_SB PHASE 3: RASTTK ({wildcards.genome}) ==="
        DOMAIN=$(awk -F'\\t' 'NR==1{{for(i=1;i<=NF;i++) if($i=="GTDBTK_domain") c=i}} NR==2{{print $c}}' {input.gtdbtk_results})
        # No silent default here on purpose: split_gtdbtk_batch_per_genome
        # (or the output_cache restore that stands in for it) always writes
        # this genome's OWN translation_table.tsv, in THIS run's own
        # output_dir, with a real header -- if "translation_table" isn't
        # found or the row is empty, that's a genuine upstream problem
        # worth failing loudly on, not a case to silently guess code 2 for
        # (margie_sb.smk's run_gtdbtk_batch comment documents exactly this
        # kind of silent-wrong-genetic-code failure mode from before the
        # ANI-only fast path was disabled).
        GCODE=$(awk -F'\t' 'NR==1{{for(i=1;i<=NF;i++){{gsub(/\r/,"",$i); if($i=="translation_table") c=i}}}} NR>1{{gsub(/\r/,"",$c); if(c && $c!=""){{print $c; exit}}}}' {input.translation_table})
        if [[ -z "$GCODE" ]]; then
            echo "ERROR: could not read a genetic code for {wildcards.genome} from {input.translation_table} -- refusing to guess" >&2
            exit 1
        fi

        # 1x, not 2x: SLURM itself already guarantees no real holder can
        # still be running past its own walltime allocation, so any lock
        # older than that is a hard-killed holder whose EXIT trap never ran.
        STALE_SECONDS=$(( {resources.runtime} * 60 ))
        while ! mkdir "{params.lock}" 2>/dev/null; do
            if [[ -f "{params.lock}/acquired_at" ]]; then
                AGE=$(( $(date +%s) - $(cat "{params.lock}/acquired_at" 2>/dev/null || echo 0) ))
                if [[ $AGE -gt $STALE_SECONDS ]]; then
                    echo "WARNING: stale RASTtk/BV-BRC lock (age ${{AGE}}s) -- assuming the holder was killed, taking over" >&2
                    rm -rf "{params.lock}"
                    continue
                fi
            fi
            sleep 10
        done
        date +%s > "{params.lock}/acquired_at"
        trap 'rm -rf "{params.lock}"' EXIT

        /usr/local/bin/run -i {input.fasta} -o {params.output_dir} -t {threads} -d {params.db} --scientific {wildcards.genome} --domain "$DOMAIN" --genetic-code "$GCODE"
        cp $(find {params.output_dir} -path "*/processed/rast*.tsv") {output.results}
        cp $(find {params.output_dir} -name "genome.faa") {output.faa}
        cp $(find {params.output_dir} -name "genome.gff") {output.gff}
        cp $(find {params.output_dir} -type d -name gene_calls -print -quit)/* $(dirname {output.results})/ || true
        echo "rasttk complete for {wildcards.genome}" > {output.tkn}
        """


rule run_prodigal:
    """Phase3: Prodigal gene calls for PRODIGAL_GENOMES (no known domain or genetic code).

    Writes RASTtk's layout (rast.tsv/.faa/.gff, shared feature ids) without
    functional descriptions or EC numbers. -g defaults to 11.
    """
    input:
        fasta=lambda wildcards: GENOMES[wildcards.genome],
        info=GENOME_INFO
    output:
        results=RASTTK_RESULTS,
        faa=RASTTK_FAA,
        gff=RASTTK_GFF,
        tkn=RASTTK_COMPUTE_TOKEN
    wildcard_constraints:
        genome=_one_of(PRODIGAL_GENOMES)
    threads: rc('prodigal.threads', 1, config=config)
    resources:
        mem_mb=rc('prodigal.mem_mb', 4000, config=config),
        runtime=runtime_min('prodigal.runtime', 60, config=config)
    params:
        output_dir=lambda wildcards: PRODIGAL_CONTAINER_OUTPUTS.format(genome=wildcards.genome)
    container: sif_path('prodigal.sif', config=config, workflow_id='margie_sb')
    shell:
        """
        echo "=== MARGIE_SB PHASE 3: PRODIGAL ({wildcards.genome}) ==="
        DOMAIN=$(awk -F'\\t' 'NR==1{{for(i=1;i<=NF;i++) if($i=="GTDBTK_domain") c=i}} NR==2{{print $c}}' {input.info})
        GCODE=$(awk -F'\t' 'NR==1{{for(i=1;i<=NF;i++){{gsub(/\r/,"",$i); if($i=="translation_table") c=i}}}} NR>1{{gsub(/\r/,"",$c); if(c && $c!=""){{print $c; exit}}}}' {input.info})
        /usr/local/bin/run -i {input.fasta} -o {params.output_dir} -g "${{GCODE:-11}}" --domain "${{DOMAIN:-Unknown}}" --organism-name {wildcards.genome} --force
        mkdir -p $(dirname {output.results})
        cp {params.output_dir}/gene_calls/* $(dirname {output.results})/
        cp {params.output_dir}/gene_calls/genome.faa {output.faa}
        cp {params.output_dir}/gene_calls/genome.gff {output.gff}
        cp {params.output_dir}/processed/prodigal_gene_calls.tsv {output.results}
        echo prodigal > $(dirname {output.results})/gene_caller.txt
        echo "prodigal complete for {wildcards.genome}" > {output.tkn}
        """


rule load_rasttk_to_db:
    """Loads RASTtk results into main_database."""
    input:
        results=RASTTK_RESULTS
    output:
        tkn=RASTTK_TOKEN
    params:
        db=MAIN_DATABASE,
        script=LOAD_SCRIPT
    shell:
        """
        {LOADER_PYTHON} {params.script} tsv {input.results} {params.db} rasttk --token {output.tkn}
        """


rule run_cog:
    """Phase4: COG functional categories with RPS-BLAST (shared phase4 entrypoint contract)."""
    input:
        faa=pc_faa('cog'),
        gtdbtk_results=GENOME_INFO
    output:
        results=PHASE4_RESULTS['cog'],
        tkn=PHASE4_COMPUTE_TOKENS['cog']
    # Priority 1 (also pfam, dbcan, geneprop) stops the scheduler's tie-break
    # from starving these rules of margie_sb_phase4_slot.
    priority: 1
    threads: rc('cog.threads', 8, config=config)
    resources:
        mem_mb=rc('cog.mem_mb', 4000, config=config),
        runtime=runtime_min('cog.runtime', 60, config=config),
        margie_sb_phase4_slot=1
    params:
        pc=pc_dir('cog'),
        merge=PC_MERGE,
        output_dir=lambda wildcards: rc('cog.output_dir', f"{CONTAINER_OUTPUTS_PREFIX}cog".format(genome=wildcards.genome), config=config),
        db=db_path('cog', config=config, workflow_id='margie_sb'),
        evalue=rc('cog.evalue', '1e-2', config=config)
    container: sif_path('cog.sif', config=config, workflow_id='margie_sb')
    shell:
        """
        echo "=== MARGIE_SB PHASE 4: COG ({wildcards.genome}) ==="
        DOMAIN=$(awk -F'\\t' 'NR==1{{for(i=1;i<=NF;i++) if($i=="GTDBTK_domain") c=i}} NR==2{{print $c}}' {input.gtdbtk_results})
        if [ -s {input.faa} ]; then
            /usr/local/bin/run -i {input.faa} -o {params.output_dir} -d {params.db} -t {threads} -e {params.evalue} --organism-name {wildcards.genome} --domain "$DOMAIN"
            SRC={params.output_dir}/{wildcards.genome}
        else
            SRC=/nonexistent/protein-cache  # every protein's rows come from the cache
        fi
        sh {params.merge} "$SRC/processed/cog_results.tsv" {output.results} {params.pc} cog_results.tsv tsv
        echo "cog complete for {wildcards.genome}" > {output.tkn}
        """


rule load_cog_to_db:
    """Loads COG results into main_database."""
    input:
        results=PHASE4_RESULTS['cog']
    output:
        tkn=PHASE4_TOKENS['cog']
    params:
        db=MAIN_DATABASE,
        script=LOAD_SCRIPT
    shell:
        """
        {LOADER_PYTHON} {params.script} tsv {input.results} {params.db} cog --token {output.tkn}
        """


rule run_pfam:
    """Phase4: Pfam domains with HMMER hmmscan --cut_ga."""
    input:
        faa=pc_faa('pfam'),
        gtdbtk_results=GENOME_INFO
    output:
        results=PHASE4_RESULTS['pfam'],
        tkn=PHASE4_COMPUTE_TOKENS['pfam']
    priority: 1  # see run_cog
    threads: rc('pfam.threads', 8, config=config)
    resources:
        mem_mb=rc('pfam.mem_mb', 8000, config=config),
        runtime=runtime_min('pfam.runtime', 60, config=config),
        margie_sb_phase4_slot=1
    params:
        pc=pc_dir('pfam'),
        merge=PC_MERGE,
        output_dir=lambda wildcards: rc('pfam.output_dir', f"{CONTAINER_OUTPUTS_PREFIX}pfam".format(genome=wildcards.genome), config=config),
        db=db_path('pfam', config=config, workflow_id='margie_sb')
    container: sif_path('pfam.sif', config=config, workflow_id='margie_sb')
    shell:
        """
        echo "=== MARGIE_SB PHASE 4: PFAM ({wildcards.genome}) ==="
        DOMAIN=$(awk -F'\\t' 'NR==1{{for(i=1;i<=NF;i++) if($i=="GTDBTK_domain") c=i}} NR==2{{print $c}}' {input.gtdbtk_results})
        if [ -s {input.faa} ]; then
            /usr/local/bin/run -i {input.faa} -o {params.output_dir} -d {params.db} -t {threads} --organism-name {wildcards.genome} --domain "$DOMAIN"
            SRC={params.output_dir}/{wildcards.genome}
        else
            SRC=/nonexistent/protein-cache  # every protein's rows come from the cache
        fi
        sh {params.merge} "$SRC/processed/pfam_results.tsv" {output.results} {params.pc} pfam_results.tsv tsv
        echo "pfam complete for {wildcards.genome}" > {output.tkn}
        """


rule load_pfam_to_db:
    """Loads Pfam results into main_database."""
    input:
        results=PHASE4_RESULTS['pfam']
    output:
        tkn=PHASE4_TOKENS['pfam']
    params:
        db=MAIN_DATABASE,
        script=LOAD_SCRIPT
    shell:
        """
        {LOADER_PYTHON} {params.script} tsv {input.results} {params.db} pfam --token {output.tkn}
        """


rule run_tigrfam:
    """Phase4: TIGRFAMs roles with HMMER hmmscan --cut_tc; also keeps the raw
    domtblout (TIGRFAM_DOMTBL) for geneprop."""
    input:
        faa=pc_faa('tigrfam'),
        gtdbtk_results=GENOME_INFO
    output:
        results=PHASE4_RESULTS['tigrfam'],
        domtbl=TIGRFAM_DOMTBL,
        tkn=PHASE4_COMPUTE_TOKENS['tigrfam']
    threads: rc('tigrfam.threads', 8, config=config)
    resources:
        mem_mb=rc('tigrfam.mem_mb', 4000, config=config),
        runtime=runtime_min('tigrfam.runtime', 60, config=config),
        margie_sb_phase4_slot=1
    params:
        pc=pc_dir('tigrfam'),
        merge=PC_MERGE,
        output_dir=lambda wildcards: rc('tigrfam.output_dir', f"{CONTAINER_OUTPUTS_PREFIX}tigrfam".format(genome=wildcards.genome), config=config),
        db=db_path('tigrfam', config=config, workflow_id='margie_sb')
    container: sif_path('tigrfam.sif', config=config, workflow_id='margie_sb')
    shell:
        """
        echo "=== MARGIE_SB PHASE 4: TIGRFAM ({wildcards.genome}) ==="
        DOMAIN=$(awk -F'\\t' 'NR==1{{for(i=1;i<=NF;i++) if($i=="GTDBTK_domain") c=i}} NR==2{{print $c}}' {input.gtdbtk_results})
        if [ -s {input.faa} ]; then
            /usr/local/bin/run -i {input.faa} -o {params.output_dir} -d {params.db} -t {threads} --organism-name {wildcards.genome} --domain "$DOMAIN"
            SRC={params.output_dir}/{wildcards.genome}
        else
            SRC=/nonexistent/protein-cache  # every protein's rows come from the cache
        fi
        sh {params.merge} "$SRC/processed/tigrfam_results.tsv" {output.results} {params.pc} tigrfam_results.tsv tsv
        sh {params.merge} "$SRC/raw/tigrfam_domtbl.out" {output.domtbl} {params.pc} tigrfam_domtbl.out domtbl
        echo "tigrfam complete for {wildcards.genome}" > {output.tkn}
        """


rule load_tigrfam_to_db:
    """Loads TIGRFAMs results into main_database."""
    input:
        results=PHASE4_RESULTS['tigrfam']
    output:
        tkn=PHASE4_TOKENS['tigrfam']
    params:
        db=MAIN_DATABASE,
        script=LOAD_SCRIPT
    shell:
        """
        {LOADER_PYTHON} {params.script} tsv {input.results} {params.db} tigrfam --token {output.tkn}
        """


rule run_merops:
    """Phase4: MEROPS peptidases with DIAMOND blastp."""
    input:
        faa=pc_faa('merops'),
        gtdbtk_results=GENOME_INFO
    output:
        results=PHASE4_RESULTS['merops'],
        tkn=PHASE4_COMPUTE_TOKENS['merops']
    threads: rc('merops.threads', 8, config=config)
    resources:
        mem_mb=rc('merops.mem_mb', 4000, config=config),
        runtime=runtime_min('merops.runtime', 30, config=config),
        margie_sb_phase4_slot=1
    params:
        pc=pc_dir('merops'),
        merge=PC_MERGE,
        output_dir=lambda wildcards: rc('merops.output_dir', f"{CONTAINER_OUTPUTS_PREFIX}merops".format(genome=wildcards.genome), config=config),
        db=db_path('merops', config=config, workflow_id='margie_sb'),
        evalue=rc('merops.evalue', '1e-5', config=config)
    container: sif_path('merops.sif', config=config, workflow_id='margie_sb')
    shell:
        """
        echo "=== MARGIE_SB PHASE 4: MEROPS ({wildcards.genome}) ==="
        DOMAIN=$(awk -F'\\t' 'NR==1{{for(i=1;i<=NF;i++) if($i=="GTDBTK_domain") c=i}} NR==2{{print $c}}' {input.gtdbtk_results})
        if [ -s {input.faa} ]; then
            /usr/local/bin/run -i {input.faa} -o {params.output_dir} -d {params.db} -t {threads} -e {params.evalue} --organism-name {wildcards.genome} --domain "$DOMAIN"
            SRC={params.output_dir}/{wildcards.genome}
        else
            SRC=/nonexistent/protein-cache  # every protein's rows come from the cache
        fi
        sh {params.merge} "$SRC/processed/merops_results.tsv" {output.results} {params.pc} merops_results.tsv tsv
        echo "merops complete for {wildcards.genome}" > {output.tkn}
        """


rule load_merops_to_db:
    """Loads MEROPS results into main_database."""
    input:
        results=PHASE4_RESULTS['merops']
    output:
        tkn=PHASE4_TOKENS['merops']
    params:
        db=MAIN_DATABASE,
        script=LOAD_SCRIPT
    shell:
        """
        {LOADER_PYTHON} {params.script} tsv {input.results} {params.db} merops --token {output.tkn}
        """


rule run_tcdb:
    """Phase4: TCDB transporter classes with DIAMOND blastp and a percent-identity cutoff (--id)."""
    input:
        faa=pc_faa('tcdb'),
        gtdbtk_results=GENOME_INFO
    output:
        results=PHASE4_RESULTS['tcdb'],
        tkn=PHASE4_COMPUTE_TOKENS['tcdb']
    threads: rc('tcdb.threads', 8, config=config)
    resources:
        mem_mb=rc('tcdb.mem_mb', 4000, config=config),
        runtime=runtime_min('tcdb.runtime', 30, config=config),
        margie_sb_phase4_slot=1
    params:
        pc=pc_dir('tcdb'),
        merge=PC_MERGE,
        output_dir=lambda wildcards: rc('tcdb.output_dir', f"{CONTAINER_OUTPUTS_PREFIX}tcdb".format(genome=wildcards.genome), config=config),
        db=db_path('tcdb', config=config, workflow_id='margie_sb'),
        evalue=rc('tcdb.evalue', '1e-5', config=config),
        pct_id=rc('tcdb.pct_id', '30', config=config)
    container: sif_path('tcdb.sif', config=config, workflow_id='margie_sb')
    shell:
        """
        echo "=== MARGIE_SB PHASE 4: TCDB ({wildcards.genome}) ==="
        DOMAIN=$(awk -F'\\t' 'NR==1{{for(i=1;i<=NF;i++) if($i=="GTDBTK_domain") c=i}} NR==2{{print $c}}' {input.gtdbtk_results})
        if [ -s {input.faa} ]; then
            /usr/local/bin/run -i {input.faa} -o {params.output_dir} -d {params.db} -t {threads} -e {params.evalue} --id {params.pct_id} --organism-name {wildcards.genome} --domain "$DOMAIN"
            SRC={params.output_dir}/{wildcards.genome}
        else
            SRC=/nonexistent/protein-cache  # every protein's rows come from the cache
        fi
        sh {params.merge} "$SRC/processed/tcdb_results.tsv" {output.results} {params.pc} tcdb_results.tsv tsv
        echo "tcdb complete for {wildcards.genome}" > {output.tkn}
        """


rule load_tcdb_to_db:
    """Loads TCDB results into main_database."""
    input:
        results=PHASE4_RESULTS['tcdb']
    output:
        tkn=PHASE4_TOKENS['tcdb']
    params:
        db=MAIN_DATABASE,
        script=LOAD_SCRIPT
    shell:
        """
        {LOADER_PYTHON} {params.script} tsv {input.results} {params.db} tcdb --token {output.tkn}
        """


rule run_uniprot:
    """Phase4: UniProt/Swiss-Prot homology with DIAMOND blastp."""
    input:
        faa=pc_faa('uniprot'),
        gtdbtk_results=GENOME_INFO
    output:
        results=PHASE4_RESULTS['uniprot'],
        tkn=PHASE4_COMPUTE_TOKENS['uniprot']
    threads: rc('uniprot.threads', 8, config=config)
    resources:
        mem_mb=rc('uniprot.mem_mb', 4000, config=config),
        runtime=runtime_min('uniprot.runtime', 30, config=config),
        margie_sb_phase4_slot=1
    params:
        pc=pc_dir('uniprot'),
        merge=PC_MERGE,
        output_dir=lambda wildcards: rc('uniprot.output_dir', f"{CONTAINER_OUTPUTS_PREFIX}uniprot".format(genome=wildcards.genome), config=config),
        db=db_path('uniprot', config=config, workflow_id='margie_sb'),
        evalue=rc('uniprot.evalue', '1e-5', config=config),
        pct_id=rc('uniprot.pct_id', '30', config=config)
    container: sif_path('uniprot.sif', config=config, workflow_id='margie_sb')
    shell:
        """
        echo "=== MARGIE_SB PHASE 4: UNIPROT ({wildcards.genome}) ==="
        DOMAIN=$(awk -F'\\t' 'NR==1{{for(i=1;i<=NF;i++) if($i=="GTDBTK_domain") c=i}} NR==2{{print $c}}' {input.gtdbtk_results})
        if [ -s {input.faa} ]; then
            /usr/local/bin/run -i {input.faa} -o {params.output_dir} -d {params.db} -t {threads} -e {params.evalue} --id {params.pct_id} --organism-name {wildcards.genome} --domain "$DOMAIN"
            SRC={params.output_dir}/{wildcards.genome}
        else
            SRC=/nonexistent/protein-cache  # every protein's rows come from the cache
        fi
        sh {params.merge} "$SRC/processed/uniprot_results.tsv" {output.results} {params.pc} uniprot_results.tsv tsv
        echo "uniprot complete for {wildcards.genome}" > {output.tkn}
        """


rule load_uniprot_to_db:
    """Loads UniProt results into main_database."""
    input:
        results=PHASE4_RESULTS['uniprot']
    output:
        tkn=PHASE4_TOKENS['uniprot']
    params:
        db=MAIN_DATABASE,
        script=LOAD_SCRIPT
    shell:
        """
        {LOADER_PYTHON} {params.script} tsv {input.results} {params.db} uniprot --token {output.tkn}
        """


rule run_kegg:
    """Phase4: KEGG Orthology with KofamScan (per-KO thresholds, no evalue flag)."""
    input:
        faa=pc_faa('kegg'),
        gtdbtk_results=GENOME_INFO
    output:
        results=PHASE4_RESULTS['kegg'],
        tkn=PHASE4_COMPUTE_TOKENS['kegg']
    threads: rc('kegg.threads', 8, config=config)
    resources:
        mem_mb=rc('kegg.mem_mb', 16000, config=config),
        runtime=runtime_min('kegg.runtime', 90, config=config),
        margie_sb_phase4_slot=1
    params:
        pc=pc_dir('kegg'),
        merge=PC_MERGE,
        output_dir=lambda wildcards: rc('kegg.output_dir', f"{CONTAINER_OUTPUTS_PREFIX}kegg".format(genome=wildcards.genome), config=config),
        db=db_path('kegg', config=config, workflow_id='margie_sb')
    container: sif_path('kegg.sif', config=config, workflow_id='margie_sb')
    shell:
        """
        echo "=== MARGIE_SB PHASE 4: KEGG ({wildcards.genome}) ==="
        DOMAIN=$(awk -F'\\t' 'NR==1{{for(i=1;i<=NF;i++) if($i=="GTDBTK_domain") c=i}} NR==2{{print $c}}' {input.gtdbtk_results})
        if [ -s {input.faa} ]; then
            /usr/local/bin/run -i {input.faa} -o {params.output_dir} -d {params.db} -t {threads} --organism-name {wildcards.genome} --domain "$DOMAIN"
            SRC={params.output_dir}/{wildcards.genome}
        else
            SRC=/nonexistent/protein-cache  # every protein's rows come from the cache
        fi
        sh {params.merge} "$SRC/processed/kegg_results.tsv" {output.results} {params.pc} kegg_results.tsv tsv
        echo "kegg complete for {wildcards.genome}" > {output.tkn}
        """


rule load_kegg_to_db:
    """Loads KEGG results into main_database."""
    input:
        results=PHASE4_RESULTS['kegg']
    output:
        tkn=PHASE4_TOKENS['kegg']
    params:
        db=MAIN_DATABASE,
        script=LOAD_SCRIPT
    shell:
        """
        {LOADER_PYTHON} {params.script} tsv {input.results} {params.db} kegg --token {output.tkn}
        """


rule run_eggnog:
    """Phase4: eggNOG-mapper orthology annotation."""
    input:
        faa=pc_faa('eggnog'),
        gtdbtk_results=GENOME_INFO
    output:
        results=PHASE4_RESULTS['eggnog'],
        tkn=PHASE4_COMPUTE_TOKENS['eggnog']
    threads: rc('eggnog.threads', 8, config=config)
    resources:
        mem_mb=rc('eggnog.mem_mb', 64000, config=config),
        runtime=runtime_min('eggnog.runtime', 90, config=config),
        margie_sb_phase4_slot=1
    params:
        pc=pc_dir('eggnog'),
        merge=PC_MERGE,
        output_dir=lambda wildcards: rc('eggnog.output_dir', f"{CONTAINER_OUTPUTS_PREFIX}eggnog".format(genome=wildcards.genome), config=config),
        db=db_path('eggnog', config=config, workflow_id='margie_sb')
    container: sif_path('eggnog.sif', config=config, workflow_id='margie_sb')
    shell:
        """
        echo "=== MARGIE_SB PHASE 4: EGGNOG ({wildcards.genome}) ==="
        DOMAIN=$(awk -F'\\t' 'NR==1{{for(i=1;i<=NF;i++) if($i=="GTDBTK_domain") c=i}} NR==2{{print $c}}' {input.gtdbtk_results})
        if [ -s {input.faa} ]; then
            /usr/local/bin/run -i {input.faa} -o {params.output_dir} -d {params.db} -t {threads} --organism-name {wildcards.genome} --domain "$DOMAIN"
            SRC={params.output_dir}/{wildcards.genome}
        else
            SRC=/nonexistent/protein-cache  # every protein's rows come from the cache
        fi
        sh {params.merge} "$SRC/processed/eggnog_results.tsv" {output.results} {params.pc} eggnog_results.tsv tsv
        echo "eggnog complete for {wildcards.genome}" > {output.tkn}
        """


rule load_eggnog_to_db:
    """Loads eggNOG results into main_database."""
    input:
        results=PHASE4_RESULTS['eggnog']
    output:
        tkn=PHASE4_TOKENS['eggnog']
    params:
        db=MAIN_DATABASE,
        script=LOAD_SCRIPT
    shell:
        """
        {LOADER_PYTHON} {params.script} tsv {input.results} {params.db} eggnog --token {output.tkn}
        """


rule run_dbcan:
    """Phase4: dbCAN CAZymes from DIAMOND, HMMER and sub-family HMM consensus."""
    input:
        faa=pc_faa('dbcan'),
        gtdbtk_results=GENOME_INFO
    output:
        results=PHASE4_RESULTS['dbcan'],
        tkn=PHASE4_COMPUTE_TOKENS['dbcan']
    priority: 1  # see run_cog
    threads: rc('dbcan.threads', 8, config=config)
    resources:
        mem_mb=rc('dbcan.mem_mb', 16000, config=config),
        runtime=runtime_min('dbcan.runtime', 60, config=config),
        margie_sb_phase4_slot=1
    params:
        pc=pc_dir('dbcan'),
        merge=PC_MERGE,
        output_dir=lambda wildcards: rc('dbcan.output_dir', f"{CONTAINER_OUTPUTS_PREFIX}dbcan".format(genome=wildcards.genome), config=config),
        db=db_path('dbcan', config=config, workflow_id='margie_sb')
    container: sif_path('dbcan.sif', config=config, workflow_id='margie_sb')
    shell:
        """
        echo "=== MARGIE_SB PHASE 4: DBCAN ({wildcards.genome}) ==="
        DOMAIN=$(awk -F'\\t' 'NR==1{{for(i=1;i<=NF;i++) if($i=="GTDBTK_domain") c=i}} NR==2{{print $c}}' {input.gtdbtk_results})
        if [ -s {input.faa} ]; then
            /usr/local/bin/run -i {input.faa} -o {params.output_dir} -d {params.db} -t {threads} --organism-name {wildcards.genome} --domain "$DOMAIN"
            SRC={params.output_dir}/{wildcards.genome}
        else
            SRC=/nonexistent/protein-cache  # every protein's rows come from the cache
        fi
        sh {params.merge} "$SRC/processed/dbcan_results.tsv" {output.results} {params.pc} dbcan_results.tsv tsv
        echo "dbcan complete for {wildcards.genome}" > {output.tkn}
        """


rule load_dbcan_to_db:
    """Loads dbCAN results into main_database."""
    input:
        results=PHASE4_RESULTS['dbcan']
    output:
        tkn=PHASE4_TOKENS['dbcan']
    params:
        db=MAIN_DATABASE,
        script=LOAD_SCRIPT
    shell:
        """
        {LOADER_PYTHON} {params.script} tsv {input.results} {params.db} dbcan --token {output.tkn}
        """


rule run_pgap:
    """Phase4: PGAP HMMs with HMMER hmmscan --cut_tc against NCBI's hmm_PGAP.LIB."""
    input:
        faa=pc_faa('pgap'),
        gtdbtk_results=GENOME_INFO
    output:
        results=PHASE4_RESULTS['pgap'],
        tkn=PHASE4_COMPUTE_TOKENS['pgap']
    threads: rc('pgap.threads', 8, config=config)
    resources:
        mem_mb=rc('pgap.mem_mb', 12000, config=config),
        runtime=runtime_min('pgap.runtime', 60, config=config),
        margie_sb_phase4_slot=1
    params:
        pc=pc_dir('pgap'),
        merge=PC_MERGE,
        output_dir=lambda wildcards: rc('pgap.output_dir', f"{CONTAINER_OUTPUTS_PREFIX}pgap".format(genome=wildcards.genome), config=config),
        db=db_path('pgap', config=config, workflow_id='margie_sb')
    container: sif_path('pgap.sif', config=config, workflow_id='margie_sb')
    shell:
        """
        echo "=== MARGIE_SB PHASE 4: PGAP ({wildcards.genome}) ==="
        DOMAIN=$(awk -F'\\t' 'NR==1{{for(i=1;i<=NF;i++) if($i=="GTDBTK_domain") c=i}} NR==2{{print $c}}' {input.gtdbtk_results})
        if [ -s {input.faa} ]; then
            /usr/local/bin/run -i {input.faa} -o {params.output_dir} -d {params.db} -t {threads} --organism-name {wildcards.genome} --domain "$DOMAIN"
            SRC={params.output_dir}/{wildcards.genome}
        else
            SRC=/nonexistent/protein-cache  # every protein's rows come from the cache
        fi
        sh {params.merge} "$SRC/processed/pgap_results.tsv" {output.results} {params.pc} pgap_results.tsv tsv
        echo "pgap complete for {wildcards.genome}" > {output.tkn}
        """


rule load_pgap_to_db:
    """Loads PGAP results into main_database."""
    input:
        results=PHASE4_RESULTS['pgap']
    output:
        tkn=PHASE4_TOKENS['pgap']
    params:
        db=MAIN_DATABASE,
        script=LOAD_SCRIPT
    shell:
        """
        {LOADER_PYTHON} {params.script} tsv {input.results} {params.db} pgap --token {output.tkn}
        """


rule run_interpro:
    """Phase4: InterProScan over the active analyses (-a), writing the unified table
    plus one per-database TSV per analysis (always written, even with no hits).

    threads and mem_mb are sized for the default four analyses on the cpu partition.
    """
    input:
        faa=pc_faa('interpro'),
        gtdbtk_results=GENOME_INFO
    output:
        results=PHASE4_RESULTS['interpro'],
        perdb=list(INTERPRO_PERDB_RESULTS.values()),
        tkn=PHASE4_COMPUTE_TOKENS['interpro']
    threads: rc('interpro.threads', 32, config=config)
    resources:
        mem_mb=rc('interpro.mem_mb', 48000, config=config),
        runtime=runtime_min('interpro.runtime', 300, config=config),
        # Default cpu partition: highmem requires at least 64 cores. The group's
        # summed mem_mb must fit one cpu node, so the load rules set small mem_mb.
        margie_sb_phase4_slot=1
    params:
        pc=pc_dir('interpro'),
        merge=PC_MERGE,
        output_dir=lambda wildcards: rc('interpro.output_dir', f"{CONTAINER_OUTPUTS_PREFIX}interpro".format(genome=wildcards.genome), config=config),
        db=db_path('interpro', config=config, workflow_id='margie_sb'),
        apps=rc('interpro.applications', INTERPRO_DEFAULT_APPS, config=config),
        db_basenames=" ".join(INTERPRO_DB_BASENAMES)
    container: sif_path('interpro.sif', config=config, workflow_id='margie_sb')
    shell:
        """
        echo "=== MARGIE_SB PHASE 4: INTERPRO ({wildcards.genome}) ==="
        DOMAIN=$(awk -F'\\t' 'NR==1{{for(i=1;i<=NF;i++) if($i=="GTDBTK_domain") c=i}} NR==2{{print $c}}' {input.gtdbtk_results})
        OUT_DIR=$(dirname {output.results})
        if [ -s {input.faa} ]; then
            /usr/local/bin/run -i {input.faa} -o {params.output_dir} -d {params.db} -t {threads} -a {params.apps} --organism-name {wildcards.genome} --domain "$DOMAIN" 2>&1 | tee "$OUT_DIR/interpro_container.log"
            SRC={params.output_dir}/{wildcards.genome}
        else
            SRC=/nonexistent/protein-cache  # every protein's rows come from the cache
        fi
        sh {params.merge} "$SRC/processed/interpro_results.tsv" {output.results} {params.pc} interpro_results.tsv tsv
        for db in {params.db_basenames}; do
            # A database with no output at all leaves an empty stub, as before.
            sh {params.merge} "$SRC/processed/interpro_${{db}}_results.tsv" "$OUT_DIR/interpro_${{db}}_results.tsv" {params.pc} "interpro_${{db}}_results.tsv" "tsv?"
        done
        echo "interpro complete for {wildcards.genome}" > {output.tkn}
        """


rule load_interpro_to_db:
    """Loads the unified InterProScan table into main_database.

    Small explicit mem_mb keeps the interpro group's summed request within one node.
    """
    input:
        results=PHASE4_RESULTS['interpro']
    output:
        tkn=PHASE4_TOKENS['interpro']
    resources:
        mem_mb=rc('interpro.load_mem_mb', 2000, config=config)
    params:
        db=MAIN_DATABASE,
        script=LOAD_SCRIPT
    shell:
        """
        {LOADER_PYTHON} {params.script} tsv {input.results} {params.db} interpro --token {output.tkn}
        """


rule load_interpro_perdb_to_db:
    """Loads one InterProScan per-database TSV into its own table (interpro_{db})."""
    input:
        results=lambda wildcards: INTERPRO_PERDB_RESULTS[wildcards.db].format(genome=wildcards.genome)
    output:
        tkn=INTERPRO_PERDB_TOKEN_PATTERN
    wildcard_constraints:
        db="|".join(INTERPRO_DB_BASENAMES)
    resources:
        mem_mb=rc('interpro.load_mem_mb', 2000, config=config)
    params:
        db_file=MAIN_DATABASE,
        script=LOAD_SCRIPT,
        table_name=lambda wildcards: f"interpro_{wildcards.db}"
    shell:
        """
        {LOADER_PYTHON} {params.script} tsv {input.results} {params.db_file} {params.table_name} --token {output.tkn}
        """


rule run_geneprop:
    """Phase4: Genome Properties assignment from tigrfam's raw domtblout against EBI's rules."""
    input:
        faa=RASTTK_FAA,
        gtdbtk_results=GENOME_INFO,
        tigrfam_domtbl=TIGRFAM_DOMTBL
    output:
        results=PHASE4_RESULTS['geneprop'],
        tkn=PHASE4_COMPUTE_TOKENS['geneprop']
    priority: 1  # see run_cog
    threads: rc('geneprop.threads', 4, config=config)
    resources:
        mem_mb=rc('geneprop.mem_mb', 2000, config=config),
        runtime=runtime_min('geneprop.runtime', 30, config=config),
        margie_sb_phase4_slot=1
    params:
        output_dir=lambda wildcards: rc('geneprop.output_dir', f"{CONTAINER_OUTPUTS_PREFIX}geneprop".format(genome=wildcards.genome), config=config),
        db=db_path('geneprop', config=config, workflow_id='margie_sb')
    container: sif_path('geneprop.sif', config=config, workflow_id='margie_sb')
    shell:
        """
        echo "=== MARGIE_SB PHASE 4: GENEPROP ({wildcards.genome}) ==="
        DOMAIN=$(awk -F'\\t' 'NR==1{{for(i=1;i<=NF;i++) if($i=="GTDBTK_domain") c=i}} NR==2{{print $c}}' {input.gtdbtk_results})
        /usr/local/bin/run -i {input.faa} -o {params.output_dir} -d {params.db} -t {threads} --tigrfam-domtbl {input.tigrfam_domtbl} --organism-name {wildcards.genome} --domain "$DOMAIN"
        cp {params.output_dir}/{wildcards.genome}/processed/geneprop_results.tsv {output.results}
        echo "geneprop complete for {wildcards.genome}" > {output.tkn}
        """


rule load_geneprop_to_db:
    """Loads Genome Properties results into main_database."""
    input:
        results=PHASE4_RESULTS['geneprop']
    output:
        tkn=PHASE4_TOKENS['geneprop']
    params:
        db=MAIN_DATABASE,
        script=LOAD_SCRIPT
    shell:
        """
        {LOADER_PYTHON} {params.script} tsv {input.results} {params.db} geneprop --token {output.tkn}
        """


rule run_operon:
    """Phase5: UniOP operon prediction from rast.faa and rast.gff (gene order and
    intergenic distance); no database."""
    input:
        faa=RASTTK_FAA,
        gff=RASTTK_GFF,
        gtdbtk_results=GENOME_INFO
    output:
        results=OPERON_RESULTS,
        tkn=OPERON_COMPUTE_TOKEN
    threads: rc('operon.threads', 4, config=config)
    resources:
        mem_mb=rc('operon.mem_mb', 2000, config=config),
        runtime=runtime_min('operon.runtime', 30, config=config)
    params:
        output_dir=lambda wildcards: rc('operon.output_dir', f"{CONTAINER_OUTPUTS_PREFIX}operon".format(genome=wildcards.genome), config=config)
    container: sif_path('operon.sif', config=config, workflow_id='margie_sb')
    shell:
        """
        echo "=== MARGIE_SB PHASE 5: OPERON ({wildcards.genome}) ==="
        DOMAIN=$(awk -F'\\t' 'NR==1{{for(i=1;i<=NF;i++) if($i=="GTDBTK_domain") c=i}} NR==2{{print $c}}' {input.gtdbtk_results})
        /usr/local/bin/run -i {input.faa} -g {input.gff} -o {params.output_dir} -t {threads} --organism-name {wildcards.genome} --domain "$DOMAIN"
        cp {params.output_dir}/{wildcards.genome}/processed/operon_results.tsv {output.results}
        echo "operon complete for {wildcards.genome}" > {output.tkn}
        """


rule load_operon_to_db:
    """Loads operon results into main_database."""
    input:
        results=OPERON_RESULTS
    output:
        tkn=OPERON_TOKEN
    params:
        db=MAIN_DATABASE,
        script=LOAD_SCRIPT
    shell:
        """
        {LOADER_PYTHON} {params.script} tsv {input.results} {params.db} operon --token {output.tkn}
        """


rule run_phobius:
    """Phase6: Phobius transmembrane/signal-peptide prediction (single-threaded, no
    database), plus the per-protein phobius_top1.tsv summary."""
    input:
        faa=pc_faa('phobius')
    output:
        results=PHOBIUS_RESULTS,
        top1=PHOBIUS_TOP1,
        tkn=PHOBIUS_COMPUTE_TOKEN
    threads: rc('phobius.threads', 4, config=config)
    resources:
        mem_mb=rc('phobius.mem_mb', 2000, config=config),
        runtime=runtime_min('phobius.runtime', 30, config=config)
    params:
        pc=pc_dir('phobius'),
        merge=PC_MERGE,
        output_dir=lambda wildcards: rc('phobius.output_dir', f"{CONTAINER_OUTPUTS_PREFIX}phobius".format(genome=wildcards.genome), config=config)
    container: sif_path('phobius.sif', config=config, workflow_id='margie_sb')
    shell:
        """
        echo "=== MARGIE_SB PHASE 6: PHOBIUS ({wildcards.genome}) ==="
        if [ -s {input.faa} ]; then
            /usr/local/bin/run -i {input.faa} -o {params.output_dir} -t {threads} --organism-name {wildcards.genome}
            SRC={params.output_dir}/{wildcards.genome}
        else
            SRC=/nonexistent/protein-cache  # every protein's rows come from the cache
        fi
        sh {params.merge} "$SRC/processed/phobius_results.tsv" {output.results} {params.pc} phobius_results.tsv tsv
        sh {params.merge} "$SRC/processed/phobius_top1.tsv" {output.top1} {params.pc} phobius_top1.tsv tsv
        echo "phobius complete for {wildcards.genome}" > {output.tkn}
        """


rule load_phobius_to_db:
    """Loads Phobius results into main_database."""
    input:
        results=PHOBIUS_RESULTS
    output:
        tkn=PHOBIUS_TOKEN
    params:
        db=MAIN_DATABASE,
        script=LOAD_SCRIPT
    shell:
        """
        {LOADER_PYTHON} {params.script} tsv {input.results} {params.db} phobius --token {output.tkn}
        """


rule run_tmbed:
    """Phase6: TMbed transmembrane prediction on CPU.

    Uses apptainer exec directly to bind db/tmbed/t5 and cnn over the package's
    hardcoded models/ paths, which ignore HF_HOME and -d.
    """
    input:
        faa=pc_faa('tmbed')
    output:
        results=TMBED_RESULTS,
        tkn=TMBED_COMPUTE_TOKEN
    threads: rc('tmbed.threads', 4, config=config)
    resources:
        mem_mb=rc('margie_sb.tmbed.mem_mb', rc('tmbed.mem_mb', 32000, config=config), config=config),
        runtime=runtime_min('margie_sb.tmbed.runtime', rc('tmbed.runtime', 240, config=config), config=config)
    params:
        pc=pc_dir('tmbed'),
        merge=PC_MERGE,
        output_dir=lambda wildcards: rc('tmbed.output_dir', f"{CONTAINER_OUTPUTS_PREFIX}tmbed".format(genome=wildcards.genome), config=config),
        model_dir=db_path('tmbed', config=config, workflow_id='margie_sb'),
        sif=sif_path('tmbed.sif', config=config, workflow_id='margie_sb')
    shell:
        """
        echo "=== MARGIE_SB PHASE 6: TMBED ({wildcards.genome}) ==="
        if [ -s {input.faa} ]; then
            apptainer exec \
                -B {params.model_dir}/t5:/usr/local/lib/python3.11/site-packages/tmbed/models/t5 \
                -B {params.model_dir}/cnn:/usr/local/lib/python3.11/site-packages/tmbed/models/cnn \
                {params.sif} \
                /usr/local/bin/run -i {input.faa} -o {params.output_dir} -t {threads} --organism-name {wildcards.genome}
            SRC={params.output_dir}/{wildcards.genome}
        else
            SRC=/nonexistent/protein-cache  # every protein's rows come from the cache
        fi
        sh {params.merge} "$SRC/processed/tmbed_results.tsv" {output.results} {params.pc} tmbed_results.tsv tsv
        echo "tmbed complete for {wildcards.genome}" > {output.tkn}
        """


rule load_tmbed_to_db:
    """Loads TMbed results into main_database."""
    input:
        results=TMBED_RESULTS
    output:
        tkn=TMBED_TOKEN
    params:
        db=MAIN_DATABASE,
        script=LOAD_SCRIPT
    shell:
        """
        {LOADER_PYTHON} {params.script} tsv {input.results} {params.db} tmbed --token {output.tkn}
        """


rule run_signalp6:
    """Phase6: SignalP 6.0 via the signalp6 HPC module with --organism other and --format none.

    Runs on CPU (the module's torch is CUDA-only, the GPU nodes are AMD);
    scaling plateaus near 8 threads.
    """
    input:
        faa=pc_faa('signalp6')
    output:
        results=SIGNALP6_RESULTS,
        tkn=SIGNALP6_COMPUTE_TOKEN
    threads: rc('signalp6.threads', 8, config=config)
    resources:
        mem_mb=rc('signalp6.mem_mb', 5000, config=config),
        runtime=runtime_min('signalp6.runtime', 15, config=config)
    params:
        pc=pc_dir('signalp6'),
        merge=PC_MERGE,
        output_dir=lambda wildcards: rc('signalp6.output_dir', f"{CONTAINER_OUTPUTS_PREFIX}signalp6".format(genome=wildcards.genome), config=config),
        process_script=SIGNALP6_PROCESS_SCRIPT
    envmodules:
        "biocontainers/default",
        "signalp6/6.0-fast"
    shell:
        """
        echo "=== MARGIE_SB PHASE 6: SIGNALP6 ({wildcards.genome}) ==="
        RAW_DIR={params.output_dir}/raw
        PROCESSED_DIR={params.output_dir}/processed
        mkdir -p "$RAW_DIR" "$PROCESSED_DIR"
        if [ -s {input.faa} ]; then
            SIGNALP6_CMD="signalp6 --fastafile {input.faa} --output_dir $RAW_DIR --organism other --mode fast --format none"
            $SIGNALP6_CMD
            {LOADER_PYTHON} {params.process_script} \
                --input "$RAW_DIR/prediction_results.txt" \
                --output "$PROCESSED_DIR/signalp6_results.tsv" \
                --organism-name {wildcards.genome} \
                --tool-used "SignalP 6.0" \
                --command-used "$SIGNALP6_CMD" \
                --database-used "SignalP6 bundled model | biocontainers/default + signalp6/6.0-fast" \
                --input-path {input.faa} \
                --output-path "$RAW_DIR"
            SRC="$PROCESSED_DIR"
        else
            SRC=/nonexistent/protein-cache  # every protein's rows come from the cache
        fi
        sh {params.merge} "$SRC/signalp6_results.tsv" {output.results} {params.pc} signalp6_results.tsv tsv
        echo "signalp6 complete for {wildcards.genome}" > {output.tkn}
        """


rule load_signalp6_to_db:
    """Loads SignalP6 results into main_database."""
    input:
        results=SIGNALP6_RESULTS
    output:
        tkn=SIGNALP6_TOKEN
    params:
        db=MAIN_DATABASE,
        script=LOAD_SCRIPT
    shell:
        """
        {LOADER_PYTHON} {params.script} tsv {input.results} {params.db} signalp6 --token {output.tkn}
        """


rule run_envelope:
    """Phase7: infers monoderm vs diderm envelope from tigrfam/pgap/pfam/uniprot marker hits."""
    input:
        tigrfam=PHASE4_RESULTS['tigrfam'],
        pgap=PHASE4_RESULTS['pgap'],
        pfam=PHASE4_RESULTS['pfam'],
        uniprot=PHASE4_RESULTS['uniprot'],
        gtdbtk_results=GENOME_INFO
    output:
        results=ENVELOPE_RESULTS,
        summary=ENVELOPE_SUMMARY,
        tkn=ENVELOPE_COMPUTE_TOKEN
    threads: rc('envelope.threads', 1, config=config)
    resources:
        mem_mb=rc('envelope.mem_mb', 1000, config=config),
        runtime=runtime_min('envelope.runtime', 15, config=config)
    params:
        input_dir=lambda wildcards: GENOME_PREFIX.format(genome=wildcards.genome).rstrip('/'),
        output_dir=lambda wildcards: rc('envelope.output_dir', f"{CONTAINER_OUTPUTS_PREFIX}envelope".format(genome=wildcards.genome), config=config)
    container: sif_path('envelope.sif', config=config, workflow_id='margie_sb')
    shell:
        """
        echo "=== MARGIE_SB PHASE 7: ENVELOPE ({wildcards.genome}) ==="
        DOMAIN=$(awk -F'\\t' 'NR==1{{for(i=1;i<=NF;i++) if($i=="GTDBTK_domain") c=i}} NR==2{{print $c}}' {input.gtdbtk_results})
        /usr/local/bin/run -i {params.input_dir} -o {params.output_dir} -t {threads} --organism-name {wildcards.genome} --domain "$DOMAIN"
        cp {params.output_dir}/processed/envelope_results.tsv {output.results}
        cp {params.output_dir}/processed/envelope_summary.tsv {output.summary}
        echo "envelope complete for {wildcards.genome}" > {output.tkn}
        """


rule load_envelope_to_db:
    """Loads envelope results into main_database."""
    input:
        results=ENVELOPE_RESULTS
    output:
        tkn=ENVELOPE_TOKEN
    params:
        db=MAIN_DATABASE,
        script=LOAD_SCRIPT
    shell:
        """
        {LOADER_PYTHON} {params.script} tsv {input.results} {params.db} envelope --token {output.tkn}
        """


rule run_deepsig:
    """Phase8: DeepSig signal peptides, with the required -k mapped from ENVELOPE_SUMMARY."""
    input:
        faa=RASTTK_FAA,
        envelope_summary=ENVELOPE_SUMMARY
    output:
        results=DEEPSIG_RESULTS,
        tkn=DEEPSIG_COMPUTE_TOKEN
    threads: rc('deepsig.threads', 4, config=config)
    resources:
        mem_mb=rc('deepsig.mem_mb', 2000, config=config),
        runtime=runtime_min('deepsig.runtime', 30, config=config),
        margie_sb_phase8_slot=1
    params:
        output_dir=lambda wildcards: rc('deepsig.output_dir', f"{CONTAINER_OUTPUTS_PREFIX}deepsig".format(genome=wildcards.genome), config=config)
    container: sif_path('deepsig.sif', config=config, workflow_id='margie_sb')
    shell:
        """
        echo "=== MARGIE_SB PHASE 8: DEEPSIG ({wildcards.genome}) ==="
        ENVTYPE=$(awk -F'\\t' 'NR==1{{for(i=1;i<=NF;i++) if($i=="envelope_type") c=i}} NR==2{{print $c}}' {input.envelope_summary})
        case "$ENVTYPE" in
            diderm*)   ORGCLASS=GRAM- ;;
            monoderm*) ORGCLASS=GRAM+ ;;
            archaea)   ORGCLASS=ARCH ;;
            *)         ORGCLASS=GRAM- ;;
        esac
        echo "[deepsig] envelope_type=$ENVTYPE -> -k $ORGCLASS"
        /usr/local/bin/run -i {input.faa} -o {params.output_dir} -t {threads} -k "$ORGCLASS" --organism-name {wildcards.genome}
        cp {params.output_dir}/{wildcards.genome}/processed/deepsig_results.tsv {output.results}
        echo "deepsig complete for {wildcards.genome}" > {output.tkn}
        """


rule load_deepsig_to_db:
    """Adds the ENVELOPE_* columns to DeepSig results in place, then loads them.

    Enrichment runs here on the host, since the container's python is incompatible.
    """
    input:
        results=DEEPSIG_RESULTS,
        envelope_summary=ENVELOPE_SUMMARY
    output:
        tkn=DEEPSIG_TOKEN
    params:
        db=MAIN_DATABASE,
        script=LOAD_SCRIPT,
        enrich_script=ENRICH_SCRIPT
    shell:
        """
        {LOADER_PYTHON} {params.enrich_script} --input {input.results} --envelope-summary {input.envelope_summary} --output {input.results}
        {LOADER_PYTHON} {params.script} tsv {input.results} {params.db} deepsig --token {output.tkn}
        """


rule run_psortb:
    """Phase8: PSORTb v3 localization, with the required -k (n|p|a) mapped from ENVELOPE_SUMMARY."""
    input:
        faa=RASTTK_FAA,
        envelope_summary=ENVELOPE_SUMMARY
    output:
        results=PSORTB_RESULTS,
        tkn=PSORTB_COMPUTE_TOKEN
    threads: rc('psortb.threads', 4, config=config)
    resources:
        mem_mb=rc('psortb.mem_mb', 2000, config=config),
        runtime=runtime_min('psortb.runtime', 30, config=config),
        margie_sb_phase8_slot=1
    params:
        output_dir=lambda wildcards: rc('psortb.output_dir', f"{CONTAINER_OUTPUTS_PREFIX}psortb".format(genome=wildcards.genome), config=config)
    container: sif_path('psortb.sif', config=config, workflow_id='margie_sb')
    shell:
        """
        echo "=== MARGIE_SB PHASE 8: PSORTB ({wildcards.genome}) ==="
        ENVTYPE=$(awk -F'\\t' 'NR==1{{for(i=1;i<=NF;i++) if($i=="envelope_type") c=i}} NR==2{{print $c}}' {input.envelope_summary})
        case "$ENVTYPE" in
            diderm*)   GRAMCLASS=n ;;
            monoderm*) GRAMCLASS=p ;;
            archaea)   GRAMCLASS=a ;;
            *)         GRAMCLASS=n ;;
        esac
        echo "[psortb] envelope_type=$ENVTYPE -> -k $GRAMCLASS"
        /usr/local/bin/run -i {input.faa} -o {params.output_dir} -t {threads} -k "$GRAMCLASS" --organism-name {wildcards.genome}
        cp {params.output_dir}/{wildcards.genome}/processed/psortb_results.tsv {output.results}
        echo "psortb complete for {wildcards.genome}" > {output.tkn}
        """


rule load_psortb_to_db:
    """Adds the ENVELOPE_* columns to PSORTb results on the host, then loads them."""
    input:
        results=PSORTB_RESULTS,
        envelope_summary=ENVELOPE_SUMMARY
    output:
        tkn=PSORTB_TOKEN
    params:
        db=MAIN_DATABASE,
        script=LOAD_SCRIPT,
        enrich_script=ENRICH_SCRIPT
    shell:
        """
        {LOADER_PYTHON} {params.enrich_script} --input {input.results} --envelope-summary {input.envelope_summary} --output {input.results}
        {LOADER_PYTHON} {params.script} tsv {input.results} {params.db} psortb --token {output.tkn}
        """


rule run_signalp4:
    """Phase8: SignalP 4.1 via the signalp4 module (command `signalp`); archaea map to gram-."""
    input:
        faa=RASTTK_FAA,
        envelope_summary=ENVELOPE_SUMMARY
    output:
        results=SIGNALP4_RESULTS,
        tkn=SIGNALP4_COMPUTE_TOKEN
    threads: rc('signalp4.threads', 2, config=config)
    resources:
        mem_mb=rc('signalp4.mem_mb', 2000, config=config),
        runtime=runtime_min('signalp4.runtime', 30, config=config),
        margie_sb_phase8_slot=1
    params:
        output_dir=lambda wildcards: rc('signalp4.output_dir', f"{CONTAINER_OUTPUTS_PREFIX}signalp4".format(genome=wildcards.genome), config=config),
        process_script=SIGNALP4_PROCESS_SCRIPT
    envmodules:
        "biocontainers/default",
        "signalp4/4.1"
    shell:
        """
        echo "=== MARGIE_SB PHASE 8: SIGNALP4 ({wildcards.genome}) ==="
        ENVTYPE=$(awk -F'\\t' 'NR==1{{for(i=1;i<=NF;i++) if($i=="envelope_type") c=i}} NR==2{{print $c}}' {input.envelope_summary})
        case "$ENVTYPE" in
            diderm*)   GRAMCLASS=gram- ;;
            monoderm*) GRAMCLASS=gram+ ;;
            *)         GRAMCLASS=gram- ;;
        esac
        echo "[signalp4] envelope_type=$ENVTYPE -> -t $GRAMCLASS"
        RAW_DIR={params.output_dir}/raw
        PROCESSED_DIR={params.output_dir}/processed
        mkdir -p "$RAW_DIR" "$PROCESSED_DIR"
        SIGNALP4_CMD="signalp -f short -t $GRAMCLASS {input.faa}"
        $SIGNALP4_CMD > "$RAW_DIR/signalp4_out.txt" 2> "$RAW_DIR/signalp4_err.txt"
        {LOADER_PYTHON} {params.process_script} \
            --input "$RAW_DIR/signalp4_out.txt" \
            --output "$PROCESSED_DIR/signalp4_results.tsv" \
            --organism-name {wildcards.genome} \
            --gram-class "$GRAMCLASS" \
            --command-used "$SIGNALP4_CMD" \
            --database-used "SignalP4 bundled model | biocontainers/default + signalp4/4.1" \
            --input-path {input.faa} \
            --output-path "$RAW_DIR"
        cp "$PROCESSED_DIR/signalp4_results.tsv" {output.results}
        echo "signalp4 complete for {wildcards.genome}" > {output.tkn}
        """


rule load_signalp4_to_db:
    """Adds the ENVELOPE_* columns to SignalP4 results on the host, then loads them."""
    input:
        results=SIGNALP4_RESULTS,
        envelope_summary=ENVELOPE_SUMMARY
    output:
        tkn=SIGNALP4_TOKEN
    params:
        db=MAIN_DATABASE,
        script=LOAD_SCRIPT,
        enrich_script=ENRICH_SCRIPT
    shell:
        """
        {LOADER_PYTHON} {params.enrich_script} --input {input.results} --envelope-summary {input.envelope_summary} --output {input.results}
        {LOADER_PYTHON} {params.script} tsv {input.results} {params.db} signalp4 --token {output.tkn}
        """


# ---- adding a tool ----
# Define <TOOL>_RESULTS/<TOOL>_TOKEN from GENOME_PREFIX, add the token to
# _phase4_8_targets_for_genome (gated by run_<tool>), and write a run_<tool>
# rule (container via sif_path, database via db_path, per-genome output_dir)
# plus a load_<tool>_to_db rule calling load_to_db.py. Rules needing a
# different SLURM partition must not share a group: with their load rule.
