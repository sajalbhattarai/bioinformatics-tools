"""
All MARGIE(SB) rules. Same framework as margie.smk: container: + shell:,
one rule per tool. Containers themselves are pipeline
agnostic (each one's /usr/local/bin/run entrypoint is fully runnable on its
own, outside any pipeline) -- the only thing this file adds is the exact
shell invocation + Snakemake-native container wiring, mirroring margie.smk.
"""
import os
import sys

# Add current directory to path to import workflow_helpers
sys.path.insert(0, os.path.dirname(workflow.snakefile))
from workflow_helpers import rc, rc_bool, fixed_path, sif_path, db_path, db_token, discover_genomes, get_workflow_prefix_for, get_container_outputs_prefix_for

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
# Appends ENVELOPE_envelope_type/inference_basis/evidence_json to a phase8
# tool's results.tsv (same value on every row -- envelope's decision is
# genome-level). Used by run_deepsig/run_psortb instead of a plain cp.
ENRICH_SCRIPT = os.path.join(WORKFLOW_DIR, "enrich_with_envelope.py")

# ─────────────────────── Path Definitions ─────────────────────── #
# Single source of truth for all file paths. Change paths here, not in rules.

# Common paths
# input_fasta can be a single genome FASTA file, OR a directory of them --
# discover_genomes() (workflow_helpers.py) returns a {genome_stem: filepath}
# map either way. GENOME_PREFIX is templated with the literal "{genome}"
# wildcard string, so every path built from it (below, and in every rule's
# output:) carries that wildcard -- Snakemake then runs each rule once per
# discovered genome, substituting the real stem in for {genome} each time.
MAIN_DATABASE = rc('main_database', config=config)

GENOMES = discover_genomes(rc('input_fasta', config=config))
if not GENOMES:
    raise ValueError(f"No genome files found for input_fasta={rc('input_fasta', config=config)!r}")

GENOME_PREFIX = get_workflow_prefix_for('{genome}', config=config)
# Sibling-of-output_dir scratch root every tool's container actually writes
# to (raw/processed/pipeline-log/whatever else its own -o produces) -- kept
# separate from GENOME_PREFIX/output_dir on purpose, see
# get_container_outputs_prefix_for()'s docstring. output_dir (GENOME_PREFIX)
# stays lean: only the final per-tool results.tsv + db token Snakemake
# actually tracks as rule outputs live there.
CONTAINER_OUTPUTS_PREFIX = get_container_outputs_prefix_for('{genome}', config=config)

_OUTPUT_ROOT = rc('output_dir', '', config=config).rstrip('/')

# Quast outputs
QUAST_RESULTS = f"{GENOME_PREFIX}quast/quast.tsv"
QUAST_TOKEN = f"{GENOME_PREFIX}quast/quast_db.tkn"

# Batch-only staging/aggregation paths for QUAST.
QUAST_BATCH_PREFIX = f"{_OUTPUT_ROOT}/original_container_outputs/quast" if _OUTPUT_ROOT else "original_container_outputs/quast"
QUAST_BATCH_STAGE_DIR = f"{QUAST_BATCH_PREFIX}/stage"
QUAST_BATCH_OUTPUT_DIR = f"{QUAST_BATCH_PREFIX}/container_outputs"
QUAST_BATCH_DONE = f"{QUAST_BATCH_PREFIX}/quast_batch.done"

# GTDB-Tk outputs
GTDBTK_RESULTS = f"{GENOME_PREFIX}gtdbtk/gtdbtk_results.tsv"
# Not loaded into the database -- this is plumbing for phase3 (RASTtk),
# which needs the genome's real NCBI genetic code, not a user-facing result.
GTDBTK_TRANSLATION_TABLE = f"{GENOME_PREFIX}gtdbtk/translation_table.tsv"
GTDBTK_TOKEN = f"{GENOME_PREFIX}gtdbtk/gtdbtk_db.tkn"

# Batch-only staging/aggregation paths for GTDB-Tk. The container is far
# more efficient when it sees the whole genome set once (shared DB/index
# warm-up) so we run one batch classify_wf, then split its combined outputs
# back into per-genome files (GTDBTK_RESULTS / GTDBTK_TRANSLATION_TABLE)
# to preserve the existing downstream rule contracts.
GTDBTK_BATCH_PREFIX = f"{_OUTPUT_ROOT}/original_container_outputs/gtdbtk" if _OUTPUT_ROOT else "original_container_outputs/gtdbtk"
GTDBTK_BATCH_STAGE_DIR = f"{GTDBTK_BATCH_PREFIX}/stage"
GTDBTK_BATCH_OUTPUT_DIR = f"{GTDBTK_BATCH_PREFIX}/container_outputs"
GTDBTK_BATCH_RESULTS = f"{GTDBTK_BATCH_PREFIX}/gtdbtk_results.tsv"
GTDBTK_BATCH_TRANSLATION_TABLE = f"{GTDBTK_BATCH_PREFIX}/gtdbtk.translation_table_summary.tsv"
GTDBTK_BATCH_DONE = f"{GTDBTK_BATCH_PREFIX}/gtdbtk_batch.done"

# RASTtk outputs. Unlike every other tool, RASTtk's real gene-caller files
# are hard dependencies for other rules (phase4's 12 tools + operon all
# need rast.faa; operon also needs rast.gff) -- not just provenance -- so
# they stay in the main output_dir alongside rast.tsv, never banished to
# container_outputs. The rest of gene_calls/ (genome.fna/.ffn/.gbk, the
# organism-prefixed duplicates, etc.) gets flattened directly into this
# same rasttk/ folder too (untracked, no separate gene_calls/ subfolder --
# rast.tsv/.faa/.gff already cover the files anything downstream actually
# depends on by name).
RASTTK_RESULTS = f"{GENOME_PREFIX}rasttk/rast.tsv"
RASTTK_TOKEN = f"{GENOME_PREFIX}rasttk/rasttk_db.tkn"
RASTTK_FAA = f"{GENOME_PREFIX}rasttk/rast.faa"
RASTTK_GFF = f"{GENOME_PREFIX}rasttk/rast.gff"

# Phase4: functional annotation (12 tools). Each takes RASTTK_FAA + GTDBTK's
# domain as input and writes <tool>_results.tsv. All 12 entrypoints share
# one contract: -i <faa> -o <output_root> -d <db> -t <threads> [extras]
# --organism-name <name> --domain <domain>, writing to
# <output_root>/<organism-name>/processed/<tool>_results.tsv -- pinning
# --organism-name to {genome} makes the path predictable, so (unlike
# quast/gtdbtk/rasttk) none of these need a find to locate their output.
# Each rule's mem_mb default tracks its own db/<tool> size on disk
# (interpro ~76G, eggnog ~48G, dbcan/kegg ~7G, pgap ~5.5G, pfam ~4.5G get
# real bumps), not the query genome size.
# Each of the 12 rules also declares margie_sb_phase4_slot=1, a named
# Snakemake resource that caps how many phase4 tools run *concurrently*,
# independent of --cores/--cpus-per-task -- without it Snakemake scheduled
# 5-6 8-thread tools at once on a 32-core job, oversubscribing real cores.
# Wired end-to-end from the frontend: workflow.py's build_executable() reads
# margie_sb.phase4.max_parallel_tools (default 4) and passes it through as
# --resources margie_sb_phase4_slot=<value>.
PHASE4_TOOLS = [
    "cog", "pfam", "tigrfam", "merops", "tcdb", "uniprot",
    "kegg", "eggnog", "dbcan", "pgap", "interpro", "geneprop",
]
PHASE4_RESULTS = {t: f"{GENOME_PREFIX}{t}/{t}_results.tsv" for t in PHASE4_TOOLS}
PHASE4_TOKENS = {t: f"{GENOME_PREFIX}{t}/{t}_db.tkn" for t in PHASE4_TOOLS}
# geneprop's --tigrfam-domtbl needs tigrfam's raw (untouched) hmmscan
# domtblout, not its own normalised processed/tigrfam_results.tsv. Still
# lives in output_dir (not container_outputs) since it's a real input: to
# run_geneprop, same reasoning as RASTTK_FAA/RASTTK_GFF.
TIGRFAM_DOMTBL = f"{GENOME_PREFIX}tigrfam/tigrfam_domtbl.out"

# Interpro per-database split outputs. PHASE4_RESULTS['interpro'] above is
# the *unified* table -- one row per domain hit with the member database
# named as a row VALUE (INTERPRO_analysis), unlike every other phase4
# tool's table (database identity baked into column names, e.g. PFAM_id).
# process_interpro_raw_results.py already writes real per-database split
# TSVs to container_outputs; they just weren't declared as Snakemake
# outputs nor loaded to the db.
# Every InterPro member-database analysis this install's interproscan.sh
# has active. Excludes Phobius/SignalP_EUK/SignalP_GRAM_NEGATIVE/
# SignalP_GRAM_POSITIVE/TMHMM (deactivated in our install; we already run
# phobius/tmbed/signalp4/6 standalone anyway) and TIGRFAM (not bundled as
# a member analysis in this InterProScan version, superseded by NCBIfam
# there). Full menu to choose from, not what actually runs -- see
# INTERPRO_ANALYSIS_TO_BASENAME below for the active subset.
INTERPRO_ALL_ANALYSES = {
    "AntiFam": "antifam", "CDD": "cdd", "Coils": "coils", "FunFam": "funfam",
    "Gene3D": "gene3d", "Hamap": "hamap", "MobiDBLite": "mobidb", "NCBIfam": "ncbifam",
    "PANTHER": "panther", "Pfam": "pfam", "PIRSF": "pirsf", "PIRSR": "pirsr",
    "PRINTS": "prints", "ProSitePatterns": "prosite_patterns", "ProSiteProfiles": "prosite_profiles",
    "SFLD": "sfld", "SMART": "smart", "SUPERFAMILY": "superfamily",
}

# Active subset, config-driven via interpro.analyses (a list of display
# names from INTERPRO_ALL_ANALYSES above, e.g. ["Hamap", "Pfam"]) -- defaults
# to 4 chosen for prokaryotic relevance, avoiding overlap with the standalone
# Pfam/TIGRFAM tools elsewhere in phase4: HAMAP (curated bacterial/archaeal
# proteomes), NCBIfam (NCBI's curated prokaryotic family HMMs), CDD (NCBI's
# conserved domain database), PIRSF (whole-protein classification, Tier 3 in
# labeling's trust hierarchy). A smaller set also keeps run_interpro's
# aggregate SLURM memory request under the cpu partition's per-node ceiling.
# The other 14, not enforced, just the default rationale: Coils/MobiDBLite
# have no real e-value; AntiFam flags spurious hits rather than annotating
# real ones; PANTHER/SMART/FunFam skew eukaryote-curated; Pfam/Gene3D/
# SUPERFAMILY are general-purpose (Pfam duplicates the standalone tool);
# PRINTS/ProSitePatterns/ProSiteProfiles/PIRSR/SFLD are lower-priority for a
# prokaryote-focused panel. Set interpro.analyses in config to override.
# INTERPRO_DB_BASENAMES/INTERPRO_PERDB_RESULTS/the per-db load rule all
# derive from whatever ends up active here.
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

# Phase5: operon prediction (UniOP). Different contract from phase4's 12
# tools -- takes RASTtk's FAA *and* GFF3 together (-i <faa> -g <gff>, gene
# order/strand matters here, not just sequence), and needs no database at
# all (no -d flag in its entrypoint), unlike every phase4 tool.
OPERON_RESULTS = f"{GENOME_PREFIX}operon/operon_results.tsv"
OPERON_TOKEN = f"{GENOME_PREFIX}operon/operon_db.tkn"

# Phase6 (per workflow_registry.py's authoritative phase numbers, not just
# file order): phobius + tmbed. Envelope-independent localization/topology
# tools -- both take a single FAA, same as operon's faa-only shape, no
# domain/gram-stain needed. signalp6 is also phase6 (registry: "HPC module,
# no envelope dependency") but has zero build-here scaffolding (no
# entrypoint, no processing script to mirror) and needs Snakemake's
# envmodules: mechanism instead of container: -- bigger, separate lift,
# deliberately not wired here yet.
PHOBIUS_RESULTS = f"{GENOME_PREFIX}phobius/phobius_results.tsv"
PHOBIUS_TOKEN = f"{GENOME_PREFIX}phobius/phobius_db.tkn"
# Per-protein summary (one row per protein) alongside the per-topology-
# segment PHOBIUS_RESULTS -- plumbing, not loaded to db, same role as
# GTDBTK_TRANSLATION_TABLE/ENVELOPE_SUMMARY.
PHOBIUS_TOP1 = f"{GENOME_PREFIX}phobius/phobius_top1.tsv"

# TMbed: deep-learning transmembrane predictor. Needs a real model
# directory (ProtT5-XL-U50 encoder, ~2.25GB) -- already cached at
# db/tmbed (confirmed populated: t5/, cnn/, the HF models--... cache dir),
# resolved the normal db_path() way. CPU-only for now (no --use-gpu) since
# this test genome is tiny; --use-gpu is there if a larger genome makes
# CPU inference too slow.
TMBED_RESULTS = f"{GENOME_PREFIX}tmbed/tmbed_results.tsv"
TMBED_TOKEN = f"{GENOME_PREFIX}tmbed/tmbed_db.tkn"

# SignalP 6.0: also phase6 (registry: "HPC module, no envelope dependency"),
# but no build-here container exists -- it's an HPC environment module
# (biocontainers/default + signalp6/6.0-fast) wrapping RCAC's own pre-built
# Apptainer image, hence envmodules: instead of container: in run_signalp6
# below. --format none is required: --format txt (the default) crashes with
# "OSError: File name too long" writing a per-protein plot file named after
# the entire FASTA header. Processing script lives at
# build-here/.../signalp6/scripts/, invoked directly by filesystem path,
# same as load_to_db.py/enrich_with_envelope.py.
SIGNALP6_RESULTS = f"{GENOME_PREFIX}signalp6/signalp6_results.tsv"
SIGNALP6_TOKEN = f"{GENOME_PREFIX}signalp6/signalp6_db.tkn"
SIGNALP6_PROCESS_SCRIPT = "/scratch/negishi/bhattar3/margie/build-here/build-containers/phase7-localization-and-topology/signalp6/scripts/process_signalp6_raw_results.py"

# Phase7: envelope type inference (monoderm vs diderm). Different shape:
# -i takes a whole directory and its entrypoint recursively searches it for
# <tool>/.../processed/*.tsv across four phase4 tools (tigrfam, pgap, pfam,
# uniprot), not a single file. -i points at the genome's Snakemake output
# root so original_container_outputs stays records-only; declaring those
# four results.tsv as input: still enforces DAG ordering even though
# envelope walks the directory itself. Also unlike phase4, -o writes
# raw/+processed/ directly with no organism-name subdirectory.
ENVELOPE_RESULTS = f"{GENOME_PREFIX}envelope/envelope_results.tsv"
ENVELOPE_TOKEN = f"{GENOME_PREFIX}envelope/envelope_db.tkn"
# Genome-level diderm/monoderm decision -- always exactly one row, even
# when envelope_results.tsv has zero marker-hit rows. Plumbing for phase8's
# envelope-dependent localization tools (psortb, deepsig, signalp4), not a
# user-facing result on its own -- same role as GTDBTK_TRANSLATION_TABLE.
ENVELOPE_SUMMARY = f"{GENOME_PREFIX}envelope/envelope_summary.tsv"

# Phase8: envelope-dependent localization (psortb, deepsig, signalp4).
# DeepSig's -k GRAM-|GRAM+|ARCH flag is a hard required argument, not
# provenance -- needs ENVELOPE_SUMMARY's real envelope_type decision, mapped
# diderm-gram-negative-like -> GRAM-, monoderm-gram-positive-like -> GRAM+,
# archaea -> ARCH. margie_sb_phase8_slot mirrors PHASE4_TOOLS' slot resource;
# workflow.py doesn't wire margie_sb.phase8.max_parallel_tools through to
# --resources yet, but workflow_registry.py already declares the config
# param, ready for whenever that's added.
DEEPSIG_RESULTS = f"{GENOME_PREFIX}deepsig/deepsig_results.tsv"
DEEPSIG_TOKEN = f"{GENOME_PREFIX}deepsig/deepsig_db.tkn"

# PSORTb v3: same envelope-dependent phase8 shape as deepsig, but -k uses
# single-letter codes (n|p|a) instead of GRAM-/GRAM+/ARCH. PSORTb itself
# often exits non-zero on warnings even when it actually succeeds -- the
# entrypoint already tolerates that internally (set +e around the perl
# call), so no extra handling needed on this side.
PSORTB_RESULTS = f"{GENOME_PREFIX}psortb/psortb_results.tsv"
PSORTB_TOKEN = f"{GENOME_PREFIX}psortb/psortb_db.tkn"

# SignalP4: also phase8 (envelope-dependent), also no build-here container
# (envmodules: biocontainers/default + signalp4/4.1 -- real command after
# module load is `signalp`, not `signalp4`). Unlike signalp6, its -t only
# supports euk/gram+/gram- -- no archaea option at all -- so archaea
# genomes get mapped to gram- in run_signalp4 below (the same conservative
# default used elsewhere when there's no real answer, not a biological
# claim). No processing script existed anywhere; wrote one at
# build-here/.../signalp4/scripts/, same non-baked-in reasoning as signalp6.
SIGNALP4_RESULTS = f"{GENOME_PREFIX}signalp4/signalp4_results.tsv"
SIGNALP4_TOKEN = f"{GENOME_PREFIX}signalp4/signalp4_db.tkn"
SIGNALP4_PROCESS_SCRIPT = "/scratch/negishi/bhattar3/margie/build-here/build-containers/phase7-localization-and-topology/signalp4/scripts/process_signalp4_raw_results.py"

# Phase9 (consolidation): modular pipeline of scripts under
# workflow_tools/consolidation/ (detect-columns.py, merge-all-columns.py,
# filter-no-stat.py). No container -- bare host-side scripts, same as
# workflow_registry.py's uses_container=False already declared.
CONSOLIDATION_SCRIPTS_DIR = os.path.join(WORKFLOW_DIR, "consolidation")
CONSOLIDATION_DETECTED_COLUMNS = f"{GENOME_PREFIX}consolidation/detected-columns.json"
CONSOLIDATION_MERGED = f"{GENOME_PREFIX}consolidation/consolidated-merged-all-columns.tsv"
CONSOLIDATION_MANIFEST = f"{GENOME_PREFIX}consolidation/manifest.tsv"
CONSOLIDATION_NO_STAT = f"{GENOME_PREFIX}consolidation/consolidated-no-stat.tsv"
CONSOLIDATION_TOKEN = f"{GENOME_PREFIX}consolidation/consolidation_db.tkn"

# Phase10 (labeling): workflow_tools/labeling/ (assign-canonical-label.py,
# add-ec-consensus.py, add-operon-info.py). assign-canonical-label.py needs
# CONSOLIDATION_MERGED specifically, not the filtered view -- its own
# docstring is explicit about needing InterPro sub-database columns and
# score/threshold columns the filtered view strips out. The other two each
# need both the labeled output AND the merged table (independent derived
# views over the same two upstream files, not chained through each other).
LABELING_SCRIPTS_DIR = os.path.join(WORKFLOW_DIR, "labeling")
LABELING_LABELED = f"{GENOME_PREFIX}labeling/labeled-genes.tsv"
LABELING_EC_CONSENSUS = f"{GENOME_PREFIX}labeling/labeled-genes-ec-consensus.tsv"
LABELING_OPERON_INFO = f"{GENOME_PREFIX}labeling/labeled-genes-operon-info.tsv"
LABELING_TOKEN = f"{GENOME_PREFIX}labeling/labeling_db.tkn"

# Sequencing contract (same "block behind failures" semantics as phase4-8):
# genome G's phase9 must not start until G's phase4-8 tools complete; phase10
# depends on phase9's merged table directly via Snakemake's own input:
# dependency. rule phase4_10_one_genome chains a per-genome token (same
# shape as PHASE4_TOKENS) through the same FIFO-queue pattern
# workflow.py's _run_pipeline_batch_sequential uses for phase4-8 -- its
# Stage 2 target is now this rule instead, one stage later.


def _phase4_8_targets_for_genome(genome):
    """Every selected phase4-8 token path for ONE genome -- the same
    rc_bool(f'run_<tool>', True, ...) gating rule all uses below, just
    scoped to a single genome instead of expand()'d across every discovered
    one. Shared by rule all (expanded over every genome, for a plain
    single-invocation run) and rule phase4_8_one_genome (one genome read
    from config['target_genome']) -- workflow.py's sequential per-organism
    orchestrator (_run_pipeline_batch_sequential) runs phase1-3 breadth-
    first across every genome via rule rasttk_all, but drives phase4-8
    through phase4_8_one_genome once per genome, one at a time, in the
    order each genome's RASTtk/GTDB-Tk actually finish -- so an organism's
    local-compute work always completes fully before the next one's starts,
    while RASTtk itself (bottlenecked on BV-BRC's remote service) keeps
    running ahead independently."""
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


def _phase9_10_targets_for_genome(genome):
    """Consolidation (phase9) and labeling (phase10) tokens for ONE
    genome, gated the same way _phase4_8_targets_for_genome gates its own
    tools -- shared by rule all and rule phase4_10_one_genome below. Must
    be defined before rule all, not after -- rule all's own input: block
    calls this at parse time, same as it calls _phase4_8_targets_for_genome
    just above."""
    targets = []
    if rc_bool('run_consolidation', True, config=config):
        targets.append(CONSOLIDATION_TOKEN.format(genome=genome))
    if rc_bool('run_labeling', True, config=config):
        targets.append(LABELING_TOKEN.format(genome=genome))
    return targets


rule all:
    # Each tool's token is only required here if rc('run_<tool>', True, ...)
    # says so -- defaults to True (tool runs) so nothing changes until the
    # frontend/API actually starts sending run_<tool>=false overrides.
    # Phase4-8 gating is centralized in _phase4_8_targets_for_genome above
    # so rule phase4_8_one_genome (the sequential per-organism
    # orchestrator's Stage 2 entry point, see its docstring) can reuse the
    # exact same per-tool selection logic for a single genome.
    input:
        (expand(QUAST_TOKEN, genome=list(GENOMES.keys())) if rc_bool('run_quast', True, config=config) else []),
        (expand(GTDBTK_TOKEN, genome=list(GENOMES.keys())) if rc_bool('run_gtdbtk', True, config=config) else []),
        (expand(RASTTK_TOKEN, genome=list(GENOMES.keys())) if rc_bool('run_rasttk', True, config=config) else []),
        [_phase4_8_targets_for_genome(genome) for genome in GENOMES.keys()],
        [_phase9_10_targets_for_genome(genome) for genome in GENOMES.keys()]


rule rasttk_all:
    """Stage 1 entry point for workflow.py's sequential per-organism
    orchestrator (_run_pipeline_batch_sequential) -- phase1-3 only, every
    genome. RASTtk is bottlenecked on BV-BRC's remote service, so it (and
    the GTDB-Tk/QUAST batches it depends on) keeps running breadth-first,
    continuously, regardless of how far phase4-8 (Stage 2, one genome at a
    time -- rule phase4_8_one_genome below) has gotten.

    GTDBTK_TOKEN is requested here (gated by run_gtdbtk, same as rule
    all's own gating below) so deselecting GTDB-Tk actually does something:
    skips load_gtdbtk_to_db, the DB load. The real classification itself
    (domain + genetic code) always still runs regardless -- run_rasttk's
    own input: block needs it unconditionally, this can't be gated away
    without breaking RASTtk."""
    input:
        expand(RASTTK_TOKEN, genome=list(GENOMES.keys())),
        (expand(GTDBTK_TOKEN, genome=list(GENOMES.keys())) if rc_bool('run_gtdbtk', True, config=config) else [])


rule phase4_8_one_genome:
    """Stage 2 entry point for the same orchestrator -- every selected
    phase4-8 token for exactly ONE genome, named via
    --config target_genome=<genome>. Scoping one Snakemake invocation to
    one genome is what makes "finish this organism's local compute before
    starting the next one's" possible at all: every tool requested in ONE
    invocation still runs in parallel against the others, still capped by
    margie_sb_phase4_slot/margie_sb_phase8_slot exactly as always -- only
    the ACROSS-genome interleaving rule all's full-batch invocation would
    otherwise allow is removed by this narrower target. Kept as its own
    rule (not folded into phase4_10_one_genome below) so anything that
    still names this target specifically keeps working unchanged."""
    input:
        _phase4_8_targets_for_genome(rc('target_genome', '', config=config))


rule phase4_10_one_genome:
    """Stage 2 entry point, extended one stage further than
    phase4_8_one_genome above: every selected phase4-8 token for one genome
    plus that genome's consolidation (phase9) and labeling (phase10)
    tokens. workflow.py's _run_pipeline_batch_sequential targets this rule
    for Stage 2 now, so a genome's local-compute work, consolidation, and
    labeling all finish before the next genome's Stage 2 starts."""
    input:
        _phase4_8_targets_for_genome(rc('target_genome', '', config=config)),
        _phase9_10_targets_for_genome(rc('target_genome', '', config=config))


rule run_consolidation:
    """Phase9 (consolidation): merges every selected phase4-8 tool's own
    results.tsv for one genome into a single one-row-per-gene table, plus
    a stripped-down readable view. Host-side scripts only, no container, run
    with {LOADER_PYTHON} same as every load_*_to_db rule.

    detect-columns.py and merge-all-columns.py are independent of each
    other; filter-no-stat.py depends on merge-all-columns.py's merged
    output. All three run sequentially in one rule -- none is expensive
    enough on a single host to need its own SLURM submission."""
    input:
        phase4_8=lambda wildcards: _phase4_8_targets_for_genome(wildcards.genome),
        gtdbtk_results=GTDBTK_RESULTS
    output:
        detected_columns=CONSOLIDATION_DETECTED_COLUMNS,
        merged=CONSOLIDATION_MERGED,
        manifest=CONSOLIDATION_MANIFEST,
        no_stat=CONSOLIDATION_NO_STAT,
        tkn=CONSOLIDATION_TOKEN
    threads: rc('consolidation.threads', 1, config=config)
    resources:
        mem_mb=rc('consolidation.mem_mb', 8000, config=config),
        runtime=rc('consolidation.runtime', 30, config=config)
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
    """Phase10 (labeling): assign-canonical-label.py decides each gene's
    canonical_label by walking the trust hierarchy over CONSOLIDATION_MERGED,
    not the no-stat view, which strips the InterPro sub-database and
    score/threshold columns it needs. add-ec-consensus.py and
    add-operon-info.py are independent derived views over the labeled
    output plus the merged table, run sequentially for the same reason as
    consolidation's scripts. Host-side only, no container."""
    input:
        merged=CONSOLIDATION_MERGED,
        consolidation_tkn=CONSOLIDATION_TOKEN
    output:
        labeled=LABELING_LABELED,
        ec_consensus=LABELING_EC_CONSENSUS,
        operon_info=LABELING_OPERON_INFO,
        tkn=LABELING_TOKEN
    threads: rc('labeling.threads', 1, config=config)
    resources:
        mem_mb=rc('labeling.mem_mb', 8000, config=config),
        runtime=rc('labeling.runtime', 30, config=config)
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
        echo "labeling complete for {wildcards.genome}" > {output.tkn}
        """


rule run_quast_batch:
    """QUAST genome quality check (margie_sb phase1), batched.

    QUAST's entrypoint already accepts a genome directory and processes each
    sample independently inside one invocation. Run it once over the whole
    input set, then split/copy per-genome processed quast.tsv files to the
    workflow's stable per-genome output paths.
    """
    input:
        list(GENOMES.values())
    output:
        done=QUAST_BATCH_DONE
    group: "quast"
    threads: rc('quast.threads', 8, config=config)
    resources:
        mem_mb=rc('quast.mem_mb', 4000, config=config),
        runtime=rc('quast.runtime', 120, config=config)
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
    """Copy batched QUAST outputs back into per-genome output paths."""
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
    """Load QUAST results into SQLite database"""
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
    """GTDB-Tk taxonomic classification (margie_sb phase2), batched.

    Run the container ONCE across the full genome set so GTDB-Tk can reuse
    its heavy DB/index warm-up in one process. Then emit the combined
    processed outputs to deterministic batch files consumed by the splitter
    rule below.
    """
    input:
        list(GENOMES.values())
    output:
        results=GTDBTK_BATCH_RESULTS,
        translation_table=GTDBTK_BATCH_TRANSLATION_TABLE,
        done=GTDBTK_BATCH_DONE
    group: "gtdbtk"
    threads: rc('gtdbtk.threads', 64, config=config)
    resources:
        mem_mb=rc('gtdbtk.mem_mb', 460000, config=config),
        runtime=rc('gtdbtk.runtime', 240, config=config),
        slurm_partition=rc('gtdbtk.partition', 'highmem', config=config)
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

        # Always force full species placement (skips GTDB-Tk's ANI-only fast
        # path) so identify/MSA runs for every genome -- ANI-only genomes
        # never get an identify-stage translation table, which silently
        # produced wrong genetic codes downstream for genomes needing a
        # non-standard table (e.g. Mycoplasmatales' code 4).
        #
        # GTDBTK_PPLACER_CPUS: the container defaults pplacer to 1 CPU
        # regardless of --cpus/-t, since its memory use on GTDB-Tk's
        # bacterial reference tree scales poorly with thread count. Set to 8
        # as a moderate starting point; check peak memory via
        # `seff <slurm_job_id>` (MaxRSS) before raising further.
        GTDBTK_PLACE_SPECIES=1 GTDBTK_PPLACER_CPUS=8 /usr/local/bin/run -i {params.stage_dir} -o {params.output_dir} -d {params.db} -t {threads} --collection-name margie_sb_batch

        cp $(find {params.output_dir} -name "gtdbtk_results.tsv") {output.results}
        cp $(find {params.output_dir} -name "gtdbtk.translation_table_summary.tsv") {output.translation_table}
        touch {output.done}
        """


rule split_gtdbtk_batch_per_genome:
    """Split batched GTDB-Tk outputs into per-genome files.

    Keeps downstream phase3+ contracts unchanged: each genome still gets its
    own gtdbtk_results.tsv and translation_table.tsv under output_dir/<genome>/.
    """
    input:
        results=GTDBTK_BATCH_RESULTS,
        translation_table=GTDBTK_BATCH_TRANSLATION_TABLE,
        done=GTDBTK_BATCH_DONE
    output:
        results=expand(GTDBTK_RESULTS, genome=list(GENOMES.keys())),
        translation_tables=expand(GTDBTK_TRANSLATION_TABLE, genome=list(GENOMES.keys()))
    run:
        import csv
        from pathlib import Path

        def _genome_from_path(path_str: str) -> str:
            p = Path(path_str)
            # .../<genome>/gtdbtk/<file>
            return p.parts[-3]

        result_targets = { _genome_from_path(p): Path(p) for p in output.results }
        translation_targets = { _genome_from_path(p): Path(p) for p in output.translation_tables }

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

        for genome, target in result_targets.items():
            row = rows_by_genome.get(genome)
            if row is None:
                raise ValueError(f"Missing GTDB-Tk row for genome '{genome}' in {input.results}")
            target.parent.mkdir(parents=True, exist_ok=True)
            with target.open('w', newline='') as out_fh:
                writer = csv.DictWriter(out_fh, fieldnames=list(row.keys()), delimiter='\t')
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

        for genome, target in translation_targets.items():
            row = rows_by_genome.get(genome)
            if row is None:
                # run_gtdbtk_batch now always forces full species placement
                # (GTDBTK_PLACE_SPECIES=1), so every genome runs identify and
                # should have a row here -- a gap means something genuinely
                # went wrong (e.g. a name mismatch), not the expected
                # ANI-only-skips-identify case this used to silently paper
                # over with a hardcoded default.
                raise ValueError(f"Missing GTDB-Tk translation table row for genome '{genome}' in {input.translation_table}")
            target.parent.mkdir(parents=True, exist_ok=True)
            with target.open('w', newline='') as out_fh:
                writer = csv.DictWriter(out_fh, fieldnames=list(row.keys()), delimiter='\t')
                writer.writeheader()
                writer.writerow(row)


rule load_gtdbtk_to_db:
    """Load GTDB-Tk classification results into SQLite database. Not
    grouped with run_gtdbtk -- the SLURM executor submits a shared group:
    as ONE job, but this rule needs only default (cpu) resources while
    run_gtdbtk needs the highmem partition, and a group can't request two
    different partitions at once."""
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


rule run_rasttk:
    """RASTtk/BV-BRC structural annotation (margie_sb phase3). Depends on
    run_gtdbtk's real outputs, not just its token -- RASTtk needs the
    genome's actual NCBI genetic code (translation_table.tsv) and domain
    (Bacteria/Archaea, from gtdbtk_results.tsv's GTDBTK_domain column) to
    annotate correctly, so this is a hard biological dependency, not an
    optional ordering hint. This means run_gtdbtk executes even if a user
    sets run_gtdbtk=false while leaving run_rasttk=true -- gtdbtk's own
    load_gtdbtk_to_db (DB load) is still skippable independently, since
    nothing downstream needs that token, only the raw result files.

    --scientific is pinned to {wildcards.genome} (entrypoint.sh only
    sanitizes spaces -- underscores in genome stems pass through
    untouched) so the final per-genome dir is predictable up to the
    domain-derived suffix entrypoint.sh appends (_bact/_arch/_unknown);
    find still locates the real processed/rast*.tsv rather than guessing
    that suffix here too.
    {genome}-wildcarded, same shape as run_quast/run_gtdbtk."""
    input:
        fasta=lambda wildcards: GENOMES[wildcards.genome],
        gtdbtk_results=GTDBTK_RESULTS,
        translation_table=GTDBTK_TRANSLATION_TABLE
    output:
        results=RASTTK_RESULTS,
        faa=RASTTK_FAA,
        gff=RASTTK_GFF
    group: "rasttk"
    threads: rc('rasttk.threads', 8, config=config)
    resources:
        mem_mb=rc('rasttk.mem_mb', 8000, config=config),
        runtime=rc('rasttk.runtime', 120, config=config),
        margie_sb_phase3_slot=1
    params:
        output_dir=lambda wildcards: rc('rasttk.output_dir', f"{CONTAINER_OUTPUTS_PREFIX}rasttk".format(genome=wildcards.genome), config=config),
        db=db_path('rasttk', config=config, workflow_id='margie_sb')
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
        GCODE=$(awk -F'\t' 'NR==1{{for(i=1;i<=NF;i++) if($i=="translation_table") c=i}} NR>1{{gsub(/\r/,"",$c); if(c && $c!=""){{print $c; exit}}}}' {input.translation_table})
        if [[ -z "$GCODE" ]]; then
            echo "ERROR: could not read a genetic code for {wildcards.genome} from {input.translation_table} -- refusing to guess" >&2
            exit 1
        fi
        /usr/local/bin/run -i {input.fasta} -o {params.output_dir} -t {threads} -d {params.db} --scientific {wildcards.genome} --domain "$DOMAIN" --genetic-code "$GCODE"
        cp $(find {params.output_dir} -path "*/processed/rast*.tsv") {output.results}
        cp $(find {params.output_dir} -name "genome.faa") {output.faa}
        cp $(find {params.output_dir} -name "genome.gff") {output.gff}
        cp $(find {params.output_dir} -type d -name gene_calls -print -quit)/* $(dirname {output.results})/ || true
        """


rule load_rasttk_to_db:
    """Load RASTtk annotation results into SQLite database"""
    input:
        results=RASTTK_RESULTS
    output:
        tkn=RASTTK_TOKEN
    group: "rasttk"
    params:
        db=MAIN_DATABASE,
        script=LOAD_SCRIPT
    shell:
        """
        {LOADER_PYTHON} {params.script} tsv {input.results} {params.db} rasttk --token {output.tkn}
        """


rule run_cog:
    """COG functional category annotation (margie_sb phase4, RPS-BLAST).
    Same shape as every other phase4 tool: --organism-name pinned to
    {wildcards.genome} makes cog's own <organism>/processed/cog_results.tsv
    path fully predictable, so unlike phase1-3 no find is needed."""
    input:
        faa=RASTTK_FAA,
        gtdbtk_results=GTDBTK_RESULTS
    output:
        results=PHASE4_RESULTS['cog']
    group: "cog"
    # Confirmed via a real run's snakemake log: cog/pfam/dbcan/geneprop got
    # zero SLURM submissions across hours of runtime despite being
    # correctly selected and planned, while every other phase4 tool sharing
    # this same margie_sb_phase4_slot resource ran repeatedly -- the
    # scheduler's tie-break was consistently passing over these four
    # whenever the shared slot was contended. priority (default 0 for
    # every rule) is checked before that tie-break, so this forces the
    # scheduler to prefer these four over their default-priority siblings.
    priority: 1
    threads: rc('cog.threads', 8, config=config)
    resources:
        mem_mb=rc('cog.mem_mb', 4000, config=config),
        runtime=rc('cog.runtime', 60, config=config),
        margie_sb_phase4_slot=1
    params:
        output_dir=lambda wildcards: rc('cog.output_dir', f"{CONTAINER_OUTPUTS_PREFIX}cog".format(genome=wildcards.genome), config=config),
        db=db_path('cog', config=config, workflow_id='margie_sb'),
        evalue=rc('cog.evalue', '1e-2', config=config)
    container: sif_path('cog.sif', config=config, workflow_id='margie_sb')
    shell:
        """
        echo "=== MARGIE_SB PHASE 4: COG ({wildcards.genome}) ==="
        DOMAIN=$(awk -F'\\t' 'NR==1{{for(i=1;i<=NF;i++) if($i=="GTDBTK_domain") c=i}} NR==2{{print $c}}' {input.gtdbtk_results})
        /usr/local/bin/run -i {input.faa} -o {params.output_dir} -d {params.db} -t {threads} -e {params.evalue} --organism-name {wildcards.genome} --domain "$DOMAIN"
        cp {params.output_dir}/{wildcards.genome}/processed/cog_results.tsv {output.results}
        """


rule load_cog_to_db:
    """Load COG annotation results into SQLite database"""
    input:
        results=PHASE4_RESULTS['cog']
    output:
        tkn=PHASE4_TOKENS['cog']
    group: "cog"
    params:
        db=MAIN_DATABASE,
        script=LOAD_SCRIPT
    shell:
        """
        {LOADER_PYTHON} {params.script} tsv {input.results} {params.db} cog --token {output.tkn}
        """


rule run_pfam:
    """Pfam domain annotation (margie_sb phase4, HMMER hmmscan --cut_ga).
    Same shape as run_cog."""
    input:
        faa=RASTTK_FAA,
        gtdbtk_results=GTDBTK_RESULTS
    output:
        results=PHASE4_RESULTS['pfam']
    group: "pfam"
    priority: 1  # see run_cog's priority comment -- same scheduler-starvation fix
    threads: rc('pfam.threads', 8, config=config)
    resources:
        mem_mb=rc('pfam.mem_mb', 8000, config=config),
        runtime=rc('pfam.runtime', 60, config=config),
        margie_sb_phase4_slot=1
    params:
        output_dir=lambda wildcards: rc('pfam.output_dir', f"{CONTAINER_OUTPUTS_PREFIX}pfam".format(genome=wildcards.genome), config=config),
        db=db_path('pfam', config=config, workflow_id='margie_sb')
    container: sif_path('pfam.sif', config=config, workflow_id='margie_sb')
    shell:
        """
        echo "=== MARGIE_SB PHASE 4: PFAM ({wildcards.genome}) ==="
        DOMAIN=$(awk -F'\\t' 'NR==1{{for(i=1;i<=NF;i++) if($i=="GTDBTK_domain") c=i}} NR==2{{print $c}}' {input.gtdbtk_results})
        /usr/local/bin/run -i {input.faa} -o {params.output_dir} -d {params.db} -t {threads} --organism-name {wildcards.genome} --domain "$DOMAIN"
        cp {params.output_dir}/{wildcards.genome}/processed/pfam_results.tsv {output.results}
        """


rule load_pfam_to_db:
    """Load Pfam annotation results into SQLite database"""
    input:
        results=PHASE4_RESULTS['pfam']
    output:
        tkn=PHASE4_TOKENS['pfam']
    group: "pfam"
    params:
        db=MAIN_DATABASE,
        script=LOAD_SCRIPT
    shell:
        """
        {LOADER_PYTHON} {params.script} tsv {input.results} {params.db} pfam --token {output.tkn}
        """


rule run_tigrfam:
    """TIGRFAMs functional role annotation (margie_sb phase4, HMMER hmmscan
    --cut_tc). Same shape as run_cog. Also exposes the raw (untouched)
    domtblout as TIGRFAM_DOMTBL -- geneprop needs that raw file, not
    tigrfam's own normalised processed/tigrfam_results.tsv."""
    input:
        faa=RASTTK_FAA,
        gtdbtk_results=GTDBTK_RESULTS
    output:
        results=PHASE4_RESULTS['tigrfam'],
        domtbl=TIGRFAM_DOMTBL
    group: "tigrfam"
    threads: rc('tigrfam.threads', 8, config=config)
    resources:
        mem_mb=rc('tigrfam.mem_mb', 4000, config=config),
        runtime=rc('tigrfam.runtime', 60, config=config),
        margie_sb_phase4_slot=1
    params:
        output_dir=lambda wildcards: rc('tigrfam.output_dir', f"{CONTAINER_OUTPUTS_PREFIX}tigrfam".format(genome=wildcards.genome), config=config),
        db=db_path('tigrfam', config=config, workflow_id='margie_sb')
    container: sif_path('tigrfam.sif', config=config, workflow_id='margie_sb')
    shell:
        """
        echo "=== MARGIE_SB PHASE 4: TIGRFAM ({wildcards.genome}) ==="
        DOMAIN=$(awk -F'\\t' 'NR==1{{for(i=1;i<=NF;i++) if($i=="GTDBTK_domain") c=i}} NR==2{{print $c}}' {input.gtdbtk_results})
        /usr/local/bin/run -i {input.faa} -o {params.output_dir} -d {params.db} -t {threads} --organism-name {wildcards.genome} --domain "$DOMAIN"
        cp {params.output_dir}/{wildcards.genome}/processed/tigrfam_results.tsv {output.results}
        cp {params.output_dir}/{wildcards.genome}/raw/tigrfam_domtbl.out {output.domtbl}
        """


rule load_tigrfam_to_db:
    """Load TIGRFAMs annotation results into SQLite database"""
    input:
        results=PHASE4_RESULTS['tigrfam']
    output:
        tkn=PHASE4_TOKENS['tigrfam']
    group: "tigrfam"
    params:
        db=MAIN_DATABASE,
        script=LOAD_SCRIPT
    shell:
        """
        {LOADER_PYTHON} {params.script} tsv {input.results} {params.db} tigrfam --token {output.tkn}
        """


rule run_merops:
    """MEROPS peptidase identification (margie_sb phase4, DIAMOND blastp).
    Same shape as run_cog."""
    input:
        faa=RASTTK_FAA,
        gtdbtk_results=GTDBTK_RESULTS
    output:
        results=PHASE4_RESULTS['merops']
    group: "merops"
    threads: rc('merops.threads', 8, config=config)
    resources:
        mem_mb=rc('merops.mem_mb', 4000, config=config),
        runtime=rc('merops.runtime', 30, config=config),
        margie_sb_phase4_slot=1
    params:
        output_dir=lambda wildcards: rc('merops.output_dir', f"{CONTAINER_OUTPUTS_PREFIX}merops".format(genome=wildcards.genome), config=config),
        db=db_path('merops', config=config, workflow_id='margie_sb'),
        evalue=rc('merops.evalue', '1e-5', config=config)
    container: sif_path('merops.sif', config=config, workflow_id='margie_sb')
    shell:
        """
        echo "=== MARGIE_SB PHASE 4: MEROPS ({wildcards.genome}) ==="
        DOMAIN=$(awk -F'\\t' 'NR==1{{for(i=1;i<=NF;i++) if($i=="GTDBTK_domain") c=i}} NR==2{{print $c}}' {input.gtdbtk_results})
        /usr/local/bin/run -i {input.faa} -o {params.output_dir} -d {params.db} -t {threads} -e {params.evalue} --organism-name {wildcards.genome} --domain "$DOMAIN"
        cp {params.output_dir}/{wildcards.genome}/processed/merops_results.tsv {output.results}
        """


rule load_merops_to_db:
    """Load MEROPS annotation results into SQLite database"""
    input:
        results=PHASE4_RESULTS['merops']
    output:
        tkn=PHASE4_TOKENS['merops']
    group: "merops"
    params:
        db=MAIN_DATABASE,
        script=LOAD_SCRIPT
    shell:
        """
        {LOADER_PYTHON} {params.script} tsv {input.results} {params.db} merops --token {output.tkn}
        """


rule run_tcdb:
    """TCDB transporter classification (margie_sb phase4, DIAMOND blastp).
    Same shape as run_cog, plus a percent-identity cutoff (--id)."""
    input:
        faa=RASTTK_FAA,
        gtdbtk_results=GTDBTK_RESULTS
    output:
        results=PHASE4_RESULTS['tcdb']
    group: "tcdb"
    threads: rc('tcdb.threads', 8, config=config)
    resources:
        mem_mb=rc('tcdb.mem_mb', 4000, config=config),
        runtime=rc('tcdb.runtime', 30, config=config),
        margie_sb_phase4_slot=1
    params:
        output_dir=lambda wildcards: rc('tcdb.output_dir', f"{CONTAINER_OUTPUTS_PREFIX}tcdb".format(genome=wildcards.genome), config=config),
        db=db_path('tcdb', config=config, workflow_id='margie_sb'),
        evalue=rc('tcdb.evalue', '1e-5', config=config),
        pct_id=rc('tcdb.pct_id', '30', config=config)
    container: sif_path('tcdb.sif', config=config, workflow_id='margie_sb')
    shell:
        """
        echo "=== MARGIE_SB PHASE 4: TCDB ({wildcards.genome}) ==="
        DOMAIN=$(awk -F'\\t' 'NR==1{{for(i=1;i<=NF;i++) if($i=="GTDBTK_domain") c=i}} NR==2{{print $c}}' {input.gtdbtk_results})
        /usr/local/bin/run -i {input.faa} -o {params.output_dir} -d {params.db} -t {threads} -e {params.evalue} --id {params.pct_id} --organism-name {wildcards.genome} --domain "$DOMAIN"
        cp {params.output_dir}/{wildcards.genome}/processed/tcdb_results.tsv {output.results}
        """


rule load_tcdb_to_db:
    """Load TCDB annotation results into SQLite database"""
    input:
        results=PHASE4_RESULTS['tcdb']
    output:
        tkn=PHASE4_TOKENS['tcdb']
    group: "tcdb"
    params:
        db=MAIN_DATABASE,
        script=LOAD_SCRIPT
    shell:
        """
        {LOADER_PYTHON} {params.script} tsv {input.results} {params.db} tcdb --token {output.tkn}
        """


rule run_uniprot:
    """UniProt/Swiss-Prot homology search (margie_sb phase4, DIAMOND
    blastp). Same shape as run_tcdb."""
    input:
        faa=RASTTK_FAA,
        gtdbtk_results=GTDBTK_RESULTS
    output:
        results=PHASE4_RESULTS['uniprot']
    group: "uniprot"
    threads: rc('uniprot.threads', 8, config=config)
    resources:
        mem_mb=rc('uniprot.mem_mb', 4000, config=config),
        runtime=rc('uniprot.runtime', 30, config=config),
        margie_sb_phase4_slot=1
    params:
        output_dir=lambda wildcards: rc('uniprot.output_dir', f"{CONTAINER_OUTPUTS_PREFIX}uniprot".format(genome=wildcards.genome), config=config),
        db=db_path('uniprot', config=config, workflow_id='margie_sb'),
        evalue=rc('uniprot.evalue', '1e-5', config=config),
        pct_id=rc('uniprot.pct_id', '30', config=config)
    container: sif_path('uniprot.sif', config=config, workflow_id='margie_sb')
    shell:
        """
        echo "=== MARGIE_SB PHASE 4: UNIPROT ({wildcards.genome}) ==="
        DOMAIN=$(awk -F'\\t' 'NR==1{{for(i=1;i<=NF;i++) if($i=="GTDBTK_domain") c=i}} NR==2{{print $c}}' {input.gtdbtk_results})
        /usr/local/bin/run -i {input.faa} -o {params.output_dir} -d {params.db} -t {threads} -e {params.evalue} --id {params.pct_id} --organism-name {wildcards.genome} --domain "$DOMAIN"
        cp {params.output_dir}/{wildcards.genome}/processed/uniprot_results.tsv {output.results}
        """


rule load_uniprot_to_db:
    """Load UniProt annotation results into SQLite database"""
    input:
        results=PHASE4_RESULTS['uniprot']
    output:
        tkn=PHASE4_TOKENS['uniprot']
    group: "uniprot"
    params:
        db=MAIN_DATABASE,
        script=LOAD_SCRIPT
    shell:
        """
        {LOADER_PYTHON} {params.script} tsv {input.results} {params.db} uniprot --token {output.tkn}
        """


rule run_kegg:
    """KEGG Orthology annotation (margie_sb phase4, KofamScan). Same shape
    as run_cog (no evalue flag -- KofamScan uses its own per-KO thresholds)."""
    input:
        faa=RASTTK_FAA,
        gtdbtk_results=GTDBTK_RESULTS
    output:
        results=PHASE4_RESULTS['kegg']
    group: "kegg"
    threads: rc('kegg.threads', 8, config=config)
    resources:
        mem_mb=rc('kegg.mem_mb', 16000, config=config),
        runtime=rc('kegg.runtime', 90, config=config),
        margie_sb_phase4_slot=1
    params:
        output_dir=lambda wildcards: rc('kegg.output_dir', f"{CONTAINER_OUTPUTS_PREFIX}kegg".format(genome=wildcards.genome), config=config),
        db=db_path('kegg', config=config, workflow_id='margie_sb')
    container: sif_path('kegg.sif', config=config, workflow_id='margie_sb')
    shell:
        """
        echo "=== MARGIE_SB PHASE 4: KEGG ({wildcards.genome}) ==="
        DOMAIN=$(awk -F'\\t' 'NR==1{{for(i=1;i<=NF;i++) if($i=="GTDBTK_domain") c=i}} NR==2{{print $c}}' {input.gtdbtk_results})
        /usr/local/bin/run -i {input.faa} -o {params.output_dir} -d {params.db} -t {threads} --organism-name {wildcards.genome} --domain "$DOMAIN"
        cp {params.output_dir}/{wildcards.genome}/processed/kegg_results.tsv {output.results}
        """


rule load_kegg_to_db:
    """Load KEGG annotation results into SQLite database"""
    input:
        results=PHASE4_RESULTS['kegg']
    output:
        tkn=PHASE4_TOKENS['kegg']
    group: "kegg"
    params:
        db=MAIN_DATABASE,
        script=LOAD_SCRIPT
    shell:
        """
        {LOADER_PYTHON} {params.script} tsv {input.results} {params.db} kegg --token {output.tkn}
        """


rule run_eggnog:
    """eggNOG-mapper orthology annotation (margie_sb phase4). Same shape
    as run_cog."""
    input:
        faa=RASTTK_FAA,
        gtdbtk_results=GTDBTK_RESULTS
    output:
        results=PHASE4_RESULTS['eggnog']
    group: "eggnog"
    threads: rc('eggnog.threads', 8, config=config)
    resources:
        mem_mb=rc('eggnog.mem_mb', 64000, config=config),
        runtime=rc('eggnog.runtime', 90, config=config),
        margie_sb_phase4_slot=1
    params:
        output_dir=lambda wildcards: rc('eggnog.output_dir', f"{CONTAINER_OUTPUTS_PREFIX}eggnog".format(genome=wildcards.genome), config=config),
        db=db_path('eggnog', config=config, workflow_id='margie_sb')
    container: sif_path('eggnog.sif', config=config, workflow_id='margie_sb')
    shell:
        """
        echo "=== MARGIE_SB PHASE 4: EGGNOG ({wildcards.genome}) ==="
        DOMAIN=$(awk -F'\\t' 'NR==1{{for(i=1;i<=NF;i++) if($i=="GTDBTK_domain") c=i}} NR==2{{print $c}}' {input.gtdbtk_results})
        /usr/local/bin/run -i {input.faa} -o {params.output_dir} -d {params.db} -t {threads} --organism-name {wildcards.genome} --domain "$DOMAIN"
        cp {params.output_dir}/{wildcards.genome}/processed/eggnog_results.tsv {output.results}
        """


rule load_eggnog_to_db:
    """Load eggNOG annotation results into SQLite database"""
    input:
        results=PHASE4_RESULTS['eggnog']
    output:
        tkn=PHASE4_TOKENS['eggnog']
    group: "eggnog"
    params:
        db=MAIN_DATABASE,
        script=LOAD_SCRIPT
    shell:
        """
        {LOADER_PYTHON} {params.script} tsv {input.results} {params.db} eggnog --token {output.tkn}
        """


rule run_dbcan:
    """dbCAN CAZyme annotation (margie_sb phase4, DIAMOND + HMMER + sub-family
    HMM consensus). Same shape as run_cog."""
    input:
        faa=RASTTK_FAA,
        gtdbtk_results=GTDBTK_RESULTS
    output:
        results=PHASE4_RESULTS['dbcan']
    group: "dbcan"
    priority: 1  # see run_cog's priority comment -- same scheduler-starvation fix
    threads: rc('dbcan.threads', 8, config=config)
    resources:
        mem_mb=rc('dbcan.mem_mb', 16000, config=config),
        runtime=rc('dbcan.runtime', 60, config=config),
        margie_sb_phase4_slot=1
    params:
        output_dir=lambda wildcards: rc('dbcan.output_dir', f"{CONTAINER_OUTPUTS_PREFIX}dbcan".format(genome=wildcards.genome), config=config),
        db=db_path('dbcan', config=config, workflow_id='margie_sb')
    container: sif_path('dbcan.sif', config=config, workflow_id='margie_sb')
    shell:
        """
        echo "=== MARGIE_SB PHASE 4: DBCAN ({wildcards.genome}) ==="
        DOMAIN=$(awk -F'\\t' 'NR==1{{for(i=1;i<=NF;i++) if($i=="GTDBTK_domain") c=i}} NR==2{{print $c}}' {input.gtdbtk_results})
        /usr/local/bin/run -i {input.faa} -o {params.output_dir} -d {params.db} -t {threads} --organism-name {wildcards.genome} --domain "$DOMAIN"
        cp {params.output_dir}/{wildcards.genome}/processed/dbcan_results.tsv {output.results}
        """


rule load_dbcan_to_db:
    """Load dbCAN annotation results into SQLite database"""
    input:
        results=PHASE4_RESULTS['dbcan']
    output:
        tkn=PHASE4_TOKENS['dbcan']
    group: "dbcan"
    params:
        db=MAIN_DATABASE,
        script=LOAD_SCRIPT
    shell:
        """
        {LOADER_PYTHON} {params.script} tsv {input.results} {params.db} dbcan --token {output.tkn}
        """


rule run_pgap:
    """PGAP HMM annotation (margie_sb phase4, HMMER hmmscan --cut_tc against
    NCBI's hmm_PGAP.LIB). Same shape as run_cog."""
    input:
        faa=RASTTK_FAA,
        gtdbtk_results=GTDBTK_RESULTS
    output:
        results=PHASE4_RESULTS['pgap']
    group: "pgap"
    threads: rc('pgap.threads', 8, config=config)
    resources:
        mem_mb=rc('pgap.mem_mb', 12000, config=config),
        runtime=rc('pgap.runtime', 60, config=config),
        margie_sb_phase4_slot=1
    params:
        output_dir=lambda wildcards: rc('pgap.output_dir', f"{CONTAINER_OUTPUTS_PREFIX}pgap".format(genome=wildcards.genome), config=config),
        db=db_path('pgap', config=config, workflow_id='margie_sb')
    container: sif_path('pgap.sif', config=config, workflow_id='margie_sb')
    shell:
        """
        echo "=== MARGIE_SB PHASE 4: PGAP ({wildcards.genome}) ==="
        DOMAIN=$(awk -F'\\t' 'NR==1{{for(i=1;i<=NF;i++) if($i=="GTDBTK_domain") c=i}} NR==2{{print $c}}' {input.gtdbtk_results})
        /usr/local/bin/run -i {input.faa} -o {params.output_dir} -d {params.db} -t {threads} --organism-name {wildcards.genome} --domain "$DOMAIN"
        cp {params.output_dir}/{wildcards.genome}/processed/pgap_results.tsv {output.results}
        """


rule load_pgap_to_db:
    """Load PGAP annotation results into SQLite database"""
    input:
        results=PHASE4_RESULTS['pgap']
    output:
        tkn=PHASE4_TOKENS['pgap']
    group: "pgap"
    params:
        db=MAIN_DATABASE,
        script=LOAD_SCRIPT
    shell:
        """
        {LOADER_PYTHON} {params.script} tsv {input.results} {params.db} pgap --token {output.tkn}
        """


rule run_interpro:
    """InterProScan domain/family/GO annotation (margie_sb phase4). -a runs
    only the 4 prokaryote-relevant analyses in INTERPRO_ANALYSIS_TO_BASENAME
    (Hamap/NCBIfam/CDD/PIRSF), deliberately narrowed from the full 18 --
    see that dict's own comment for the reasoning. Produces the unified
    results.tsv plus one per-database split TSV per configured analysis --
    see INTERPRO_PERDB_RESULTS above for why those exist; the container
    always writes one file per analysis passed via -a regardless of hit
    count, so every declared per-db output is guaranteed to exist even when
    an analysis finds nothing on a given genome.
    threads=32 and mem_mb are sized for just these 4 (lighter than the
    full 18) on the cpu partition -- mem_mb is an estimate pending a real
    run's actual usage, not a measured figure; adjust interpro.mem_mb in
    config if it runs short."""
    input:
        faa=RASTTK_FAA,
        gtdbtk_results=GTDBTK_RESULTS
    output:
        results=PHASE4_RESULTS['interpro'],
        perdb=list(INTERPRO_PERDB_RESULTS.values())
    group: "interpro"
    threads: rc('interpro.threads', 32, config=config)
    resources:
        mem_mb=rc('interpro.mem_mb', 48000, config=config),
        runtime=rc('interpro.runtime', 300, config=config),
        # NOT highmem -- tried that, reverted. highmem is gtdbtk's partition
        # (above) because gtdbtk genuinely asks for 64 threads, clearing this
        # cluster's real policy: highmem is reserved for jobs needing more
        # memory than a standard node, and since memory there is allocated
        # proportional to CPU count, you must request at least 64 cores.
        # `sbatch -p highmem --cpus-per-task=19 ...` was rejected outright
        # with exactly that message -- Snakemake's log never surfaces it
        # (just "Error in group interpro"), only visible with --verbose.
        #
        # Separately, even on 'cpu' (257GB/node, per sinfo), the GROUP's
        # aggregate request used to exceed a single node's memory:
        # load_interpro_to_db and the (then 18) load_interpro_perdb_to_db
        # group-mates had no mem_mb of their own, inheriting the global
        # 16000 default, and 96000 + 19*16000 = 400000 MB blew past 257400
        # (sacct showed that group job FAILED with Elapsed=00:00:00/
        # Start=None -- SLURM rejecting an unsatisfiable allocation, not a
        # runtime OOM kill). Fixed by giving those load rules their own
        # small explicit mem_mb below; narrowing -a to 4 analyses also cut
        # run_interpro's own footprint.
        margie_sb_phase4_slot=1
    params:
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
        /usr/local/bin/run -i {input.faa} -o {params.output_dir} -d {params.db} -t {threads} -a {params.apps} --organism-name {wildcards.genome} --domain "$DOMAIN" 2>&1 | tee "$OUT_DIR/interpro_container.log"
        cp {params.output_dir}/{wildcards.genome}/processed/interpro_results.tsv {output.results}
        for db in {params.db_basenames}; do
            SRC="{params.output_dir}/{wildcards.genome}/processed/interpro_${{db}}_results.tsv"
            if [[ -f "$SRC" ]]; then
                cp "$SRC" "$OUT_DIR/interpro_${{db}}_results.tsv"
            else
                echo "[interpro] WARNING: $db produced no output file at all, writing empty stub" >&2
                touch "$OUT_DIR/interpro_${{db}}_results.tsv"
            fi
        done
        """


rule load_interpro_to_db:
    """Load InterProScan annotation results into SQLite database.
    Explicit small mem_mb (a trivial TSV->SQLite load needs nowhere near
    the 16000 global default) -- otherwise this rule's share of
    run_interpro's group resource sum adds up fast across this rule plus
    every load_interpro_perdb_to_db group-mate; see run_interpro's own
    resources comment for the incident this caused."""
    input:
        results=PHASE4_RESULTS['interpro']
    output:
        tkn=PHASE4_TOKENS['interpro']
    group: "interpro"
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
    """Load one InterProScan per-database split TSV (e.g.
    interpro_hamap_results.tsv) into its own SQLite table (e.g.
    interpro_hamap) -- one rule, matched against any of
    INTERPRO_DB_BASENAMES via the {db} wildcard, instead of one
    near-identical rule block per analysis. Same explicit small mem_mb
    as load_interpro_to_db, same reasoning."""
    input:
        results=lambda wildcards: INTERPRO_PERDB_RESULTS[wildcards.db].format(genome=wildcards.genome)
    output:
        tkn=INTERPRO_PERDB_TOKEN_PATTERN
    wildcard_constraints:
        db="|".join(INTERPRO_DB_BASENAMES)
    group: "interpro"
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
    """Genome Properties whole-genome property assignment (margie_sb
    phase4). Depends on run_tigrfam's real outputs (not just its token) --
    needs tigrfam's raw domtblout, a hard biological dependency like
    run_rasttk's dependency on gtdbtk. Lighter-weight than the HMM/BLAST
    tools (pure post-processing against EBI's genome-properties rules)."""
    input:
        faa=RASTTK_FAA,
        gtdbtk_results=GTDBTK_RESULTS,
        tigrfam_domtbl=TIGRFAM_DOMTBL
    output:
        results=PHASE4_RESULTS['geneprop']
    group: "geneprop"
    priority: 1  # see run_cog's priority comment -- same scheduler-starvation fix
    threads: rc('geneprop.threads', 4, config=config)
    resources:
        mem_mb=rc('geneprop.mem_mb', 2000, config=config),
        runtime=rc('geneprop.runtime', 30, config=config),
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
        """


rule load_geneprop_to_db:
    """Load Genome Properties results into SQLite database"""
    input:
        results=PHASE4_RESULTS['geneprop']
    output:
        tkn=PHASE4_TOKENS['geneprop']
    group: "geneprop"
    params:
        db=MAIN_DATABASE,
        script=LOAD_SCRIPT
    shell:
        """
        {LOADER_PYTHON} {params.script} tsv {input.results} {params.db} geneprop --token {output.tkn}
        """


rule run_operon:
    """Operon prediction via UniOP (margie_sb phase5). Different shape from
    every phase4 tool: needs RASTtk's GFF3 alongside its FAA (-g, gene
    order/strand matters for operon calls, not just sequence), and takes no
    database at all (no -d/db_path() -- UniOP is a pure intergenic-distance
    probabilistic model, nothing to look up). --organism-name pinned to
    {wildcards.genome} same as phase4, for the same predictable-output-path
    reason."""
    input:
        faa=RASTTK_FAA,
        gff=RASTTK_GFF,
        gtdbtk_results=GTDBTK_RESULTS
    output:
        results=OPERON_RESULTS
    group: "operon"
    threads: rc('operon.threads', 4, config=config)
    resources:
        mem_mb=rc('operon.mem_mb', 2000, config=config),
        runtime=rc('operon.runtime', 30, config=config)
    params:
        output_dir=lambda wildcards: rc('operon.output_dir', f"{CONTAINER_OUTPUTS_PREFIX}operon".format(genome=wildcards.genome), config=config)
    container: sif_path('operon.sif', config=config, workflow_id='margie_sb')
    shell:
        """
        echo "=== MARGIE_SB PHASE 5: OPERON ({wildcards.genome}) ==="
        DOMAIN=$(awk -F'\\t' 'NR==1{{for(i=1;i<=NF;i++) if($i=="GTDBTK_domain") c=i}} NR==2{{print $c}}' {input.gtdbtk_results})
        /usr/local/bin/run -i {input.faa} -g {input.gff} -o {params.output_dir} -t {threads} --organism-name {wildcards.genome} --domain "$DOMAIN"
        cp {params.output_dir}/{wildcards.genome}/processed/operon_results.tsv {output.results}
        """


rule load_operon_to_db:
    """Load operon prediction results into SQLite database"""
    input:
        results=OPERON_RESULTS
    output:
        tkn=OPERON_TOKEN
    group: "operon"
    params:
        db=MAIN_DATABASE,
        script=LOAD_SCRIPT
    shell:
        """
        {LOADER_PYTHON} {params.script} tsv {input.results} {params.db} operon --token {output.tkn}
        """


rule run_phobius:
    """Phobius combined transmembrane/signal-peptide prediction (margie_sb
    phase6). Envelope-independent, same faa-only shape as run_operon's -i
    (minus -g): single-threaded regardless of -t (entrypoint's own note),
    no database. Also exposes phobius_top1.tsv (per-protein summary)
    alongside the per-segment phobius_results.tsv."""
    input:
        faa=RASTTK_FAA
    output:
        results=PHOBIUS_RESULTS,
        top1=PHOBIUS_TOP1
    group: "phobius"
    threads: rc('phobius.threads', 4, config=config)
    resources:
        mem_mb=rc('phobius.mem_mb', 2000, config=config),
        runtime=rc('phobius.runtime', 30, config=config)
    params:
        output_dir=lambda wildcards: rc('phobius.output_dir', f"{CONTAINER_OUTPUTS_PREFIX}phobius".format(genome=wildcards.genome), config=config)
    container: sif_path('phobius.sif', config=config, workflow_id='margie_sb')
    shell:
        """
        echo "=== MARGIE_SB PHASE 6: PHOBIUS ({wildcards.genome}) ==="
        /usr/local/bin/run -i {input.faa} -o {params.output_dir} -t {threads} --organism-name {wildcards.genome}
        cp {params.output_dir}/{wildcards.genome}/processed/phobius_results.tsv {output.results}
        cp {params.output_dir}/{wildcards.genome}/processed/phobius_top1.tsv {output.top1}
        """


rule load_phobius_to_db:
    """Load Phobius prediction results into SQLite database"""
    input:
        results=PHOBIUS_RESULTS
    output:
        tkn=PHOBIUS_TOKEN
    group: "phobius"
    params:
        db=MAIN_DATABASE,
        script=LOAD_SCRIPT
    shell:
        """
        {LOADER_PYTHON} {params.script} tsv {input.results} {params.db} phobius --token {output.tkn}
        """


rule run_tmbed:
    """TMbed deep-learning transmembrane prediction (margie_sb phase6).
    Envelope-independent, same faa-only shape as run_phobius. Deliberately
    NOT using container: + the entrypoint's documented -d/HF_HOME approach
    -- confirmed by reading the actual installed code that it's wrong for
    this tmbed version: tmbed.py's load_encoder()/load_models() hardcode
    Path(__file__).parent/'models/t5' and .../'models/cnn' (the package's
    own install dir inside the image), never consulting HF_HOME, -d, or
    any CLI flag. Since the image's own filesystem is read-only, the only
    way to get the already-cached weights (db/tmbed/t5, db/tmbed/cnn --
    confirmed present: config.json, spiece.model, model.safetensors,
    cv_0-4.pt) seen at those exact paths is to bind them there directly,
    which needs a manual apptainer exec (Snakemake's container: has no
    per-rule host:container bind-path remapping). Verified directly: with
    these binds, config.json and all 5 .pt files resolve exactly where
    tmbed's hardcoded loader looks for them.
    CPU-only for now (no --use-gpu); the entrypoint supports it if a
    larger genome ever makes CPU inference too slow."""
    input:
        faa=RASTTK_FAA
    output:
        results=TMBED_RESULTS
    group: "tmbed"
    threads: rc('tmbed.threads', 4, config=config)
    resources:
        mem_mb=rc('tmbed.mem_mb', 8000, config=config),
        runtime=rc('tmbed.runtime', 60, config=config)
    params:
        output_dir=lambda wildcards: rc('tmbed.output_dir', f"{CONTAINER_OUTPUTS_PREFIX}tmbed".format(genome=wildcards.genome), config=config),
        model_dir=db_path('tmbed', config=config, workflow_id='margie_sb'),
        sif=sif_path('tmbed.sif', config=config, workflow_id='margie_sb')
    shell:
        """
        echo "=== MARGIE_SB PHASE 6: TMBED ({wildcards.genome}) ==="
        apptainer exec \
            -B {params.model_dir}/t5:/usr/local/lib/python3.11/site-packages/tmbed/models/t5 \
            -B {params.model_dir}/cnn:/usr/local/lib/python3.11/site-packages/tmbed/models/cnn \
            {params.sif} \
            /usr/local/bin/run -i {input.faa} -o {params.output_dir} -t {threads} --organism-name {wildcards.genome}
        cp {params.output_dir}/{wildcards.genome}/processed/tmbed_results.tsv {output.results}
        """


rule load_tmbed_to_db:
    """Load TMbed prediction results into SQLite database"""
    input:
        results=TMBED_RESULTS
    output:
        tkn=TMBED_TOKEN
    group: "tmbed"
    params:
        db=MAIN_DATABASE,
        script=LOAD_SCRIPT
    shell:
        """
        {LOADER_PYTHON} {params.script} tsv {input.results} {params.db} tmbed --token {output.tkn}
        """


rule run_signalp6:
    """SignalP 6.0 signal peptide prediction (margie_sb phase6).
    Envelope-independent, same faa-only shape as run_phobius/run_tmbed --
    but no container: here at all (envmodules: loads RCAC's own
    biocontainers/default + signalp6/6.0-fast HPC modules instead, which
    wrap a pre-built Apptainer image we don't own/build). --organism other
    (bacteria/archaea, not eukarya) and --format none (see
    SIGNALP6_RESULTS' comment for why -- avoids a real crash) are both
    required, not defaults to leave alone.

    GPU is not a usable speedup path here: Negishi's GPU partition is AMD
    (apptainer reports "Could not find any nv files on this host" and the
    biocontainers wrapper falls back to "Enabling AMD GPU support"), but
    this container's torch build is CUDA-only (1.11.0+cu102, zero ROCm
    support) -- confirmed via a real GPU-node test, torch.cuda.is_available()
    is False there regardless. What actually speeds this up is plain CPU
    core count: 4 threads took 12m13s real wall-clock on this genome's 530
    proteins; 10 threads took 2m49s-3m22s in two separate real runs (one on
    the GPU partition incidentally, one on the plain CPU partition -- same
    timing either way, confirming it's core count, not node type). Scaling
    plateaus around 8 effective cores (user-time/real-time ratio was ~7.7x
    at 10 allocated), hence 8 here rather than 10."""
    input:
        faa=RASTTK_FAA
    output:
        results=SIGNALP6_RESULTS
    group: "signalp6"
    threads: rc('signalp6.threads', 8, config=config)
    resources:
        mem_mb=rc('signalp6.mem_mb', 5000, config=config),
        runtime=rc('signalp6.runtime', 15, config=config)
    params:
        output_dir=lambda wildcards: rc('signalp6.output_dir', f"{CONTAINER_OUTPUTS_PREFIX}signalp6".format(genome=wildcards.genome), config=config),
        process_script=SIGNALP6_PROCESS_SCRIPT
    envmodules:
        "biocontainers/default",
        "signalp6/6.0-fast"
    shell:
        """
        echo "=== MARGIE_SB PHASE 6: SIGNALP6 ({wildcards.genome}) ==="
        mkdir -p {params.output_dir}
        SIGNALP6_CMD="signalp6 --fastafile {input.faa} --output_dir {params.output_dir} --organism other --mode fast --format none"
        $SIGNALP6_CMD
        {LOADER_PYTHON} {params.process_script} \
            --input {params.output_dir}/prediction_results.txt \
            --output {output.results} \
            --organism-name {wildcards.genome} \
            --tool-used "SignalP 6.0" \
            --command-used "$SIGNALP6_CMD" \
            --database-used "SignalP6 bundled model | biocontainers/default + signalp6/6.0-fast" \
            --input-path {input.faa} \
            --output-path {params.output_dir}
        """


rule load_signalp6_to_db:
    """Load SignalP6 prediction results into SQLite database"""
    input:
        results=SIGNALP6_RESULTS
    output:
        tkn=SIGNALP6_TOKEN
    group: "signalp6"
    params:
        db=MAIN_DATABASE,
        script=LOAD_SCRIPT
    shell:
        """
        {LOADER_PYTHON} {params.script} tsv {input.results} {params.db} signalp6 --token {output.tkn}
        """


rule run_envelope:
    """Envelope type inference (margie_sb phase6): monoderm vs diderm,
    weighted from marker hits across tigrfam/pgap/pfam/uniprot's real
    processed/ output (see ENVELOPE_RESULTS/ENVELOPE_SUMMARY's comments
    for why -i points at the genome output root rather than a single
    file). -d is accepted but silently discarded by the entrypoint
    (vestigial templating, not a real database lookup) so it's omitted
    here, same as run_operon."""
    input:
        tigrfam=PHASE4_RESULTS['tigrfam'],
        pgap=PHASE4_RESULTS['pgap'],
        pfam=PHASE4_RESULTS['pfam'],
        uniprot=PHASE4_RESULTS['uniprot'],
        gtdbtk_results=GTDBTK_RESULTS
    output:
        results=ENVELOPE_RESULTS,
        summary=ENVELOPE_SUMMARY
    group: "envelope"
    threads: rc('envelope.threads', 1, config=config)
    resources:
        mem_mb=rc('envelope.mem_mb', 1000, config=config),
        runtime=rc('envelope.runtime', 15, config=config)
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
        """


rule load_envelope_to_db:
    """Load envelope classification results into SQLite database"""
    input:
        results=ENVELOPE_RESULTS
    output:
        tkn=ENVELOPE_TOKEN
    group: "envelope"
    params:
        db=MAIN_DATABASE,
        script=LOAD_SCRIPT
    shell:
        """
        {LOADER_PYTHON} {params.script} tsv {input.results} {params.db} envelope --token {output.tkn}
        """


rule run_deepsig:
    """DeepSig signal peptide prediction (margie_sb phase8). -k is a
    required argument (entrypoint exits 1 without it) -- mapped from
    ENVELOPE_SUMMARY's real envelope_type decision, not user-supplied.
    -d is accepted by the entrypoint but silently discarded (no external
    DB needed), same as run_operon/run_phobius."""
    input:
        faa=RASTTK_FAA,
        envelope_summary=ENVELOPE_SUMMARY
    output:
        results=DEEPSIG_RESULTS
    group: "deepsig"
    threads: rc('deepsig.threads', 4, config=config)
    resources:
        mem_mb=rc('deepsig.mem_mb', 2000, config=config),
        runtime=rc('deepsig.runtime', 30, config=config),
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
        """


rule load_deepsig_to_db:
    """Enrich + load DeepSig prediction results into SQLite database.
    Enrichment (pulling in ENVELOPE_envelope_type/inference_basis/
    evidence_json) happens here, not in run_deepsig, deliberately --
    run_deepsig has container: set, so its whole shell: block executes
    via `apptainer exec ... bash -c`, where bare `python` resolves to
    deepsig.sif's own (incompatible, Python 2-era) interpreter. This rule
    has no container: (same as every other load_*_to_db rule), so it runs
    on the host under the activated venv, where `python` is real Python 3.
    Enriches DEEPSIG_RESULTS in place (overwrites with the enriched
    version) before loading -- the file in output/ ends up being the
    final, complete one either way."""
    input:
        results=DEEPSIG_RESULTS,
        envelope_summary=ENVELOPE_SUMMARY
    output:
        tkn=DEEPSIG_TOKEN
    group: "deepsig"
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
    """PSORTb v3 subcellular localization (margie_sb phase8). Same
    envelope-dependent shape as run_deepsig -- -k is required (n|p|a),
    mapped from ENVELOPE_SUMMARY's real envelope_type decision."""
    input:
        faa=RASTTK_FAA,
        envelope_summary=ENVELOPE_SUMMARY
    output:
        results=PSORTB_RESULTS
    group: "psortb"
    threads: rc('psortb.threads', 4, config=config)
    resources:
        mem_mb=rc('psortb.mem_mb', 2000, config=config),
        runtime=rc('psortb.runtime', 30, config=config),
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
        """


rule load_psortb_to_db:
    """Enrich + load PSORTb localization results into SQLite database.
    Same reasoning as load_deepsig_to_db -- enrichment must happen here
    (host-side, no container:), not in run_psortb (container: set, bare
    `python` there resolves to psortb.sif's own incompatible interpreter)."""
    input:
        results=PSORTB_RESULTS,
        envelope_summary=ENVELOPE_SUMMARY
    output:
        tkn=PSORTB_TOKEN
    group: "psortb"
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
    """SignalP 4.1 signal peptide prediction (margie_sb phase8,
    envelope-dependent). Real command after module load is `signalp`, not
    `signalp4`. -t only supports euk/gram+/gram- -- archaea genomes map to
    gram- here too (no real archaea option exists; same conservative
    default used elsewhere when there's no good answer)."""
    input:
        faa=RASTTK_FAA,
        envelope_summary=ENVELOPE_SUMMARY
    output:
        results=SIGNALP4_RESULTS
    group: "signalp4"
    threads: rc('signalp4.threads', 2, config=config)
    resources:
        mem_mb=rc('signalp4.mem_mb', 2000, config=config),
        runtime=rc('signalp4.runtime', 30, config=config),
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
        mkdir -p {params.output_dir}
        SIGNALP4_CMD="signalp -f short -t $GRAMCLASS {input.faa}"
        $SIGNALP4_CMD > {params.output_dir}/signalp4_out.txt 2> {params.output_dir}/signalp4_err.txt
        {LOADER_PYTHON} {params.process_script} \
            --input {params.output_dir}/signalp4_out.txt \
            --output {output.results} \
            --organism-name {wildcards.genome} \
            --gram-class "$GRAMCLASS" \
            --command-used "$SIGNALP4_CMD" \
            --database-used "SignalP4 bundled model | biocontainers/default + signalp4/4.1" \
            --input-path {input.faa} \
            --output-path {params.output_dir}
        """


rule load_signalp4_to_db:
    """Enrich + load SignalP4 prediction results into SQLite database.
    Same enrichment-placement reasoning as load_deepsig_to_db/
    load_psortb_to_db -- kept consistent regardless of whether the
    upstream run_* rule uses container: or envmodules:."""
    input:
        results=SIGNALP4_RESULTS,
        envelope_summary=ENVELOPE_SUMMARY
    output:
        tkn=SIGNALP4_TOKEN
    group: "signalp4"
    params:
        db=MAIN_DATABASE,
        script=LOAD_SCRIPT,
        enrich_script=ENRICH_SCRIPT
    shell:
        """
        {LOADER_PYTHON} {params.enrich_script} --input {input.results} --envelope-summary {input.envelope_summary} --output {input.results}
        {LOADER_PYTHON} {params.script} tsv {input.results} {params.db} signalp4 --token {output.tkn}
        """


# Phase9 (consolidation) rules are pending the modular pipeline rewrite --
# see the comment near the path-definitions section above.

# ═════════════════════════════════════════════════════════════════════════════
#                           NEW RULE TEMPLATE
# ═════════════════════════════════════════════════════════════════════════════
#
# Quick guide for adding new margie_sb tools, all {genome}-wildcarded so
# one Snakemake invocation processes a single genome OR a whole folder of
# them (see GENOMES/GENOME_PREFIX up top).
#
# STEP 1: Add path definitions at top, using GENOME_PREFIX (which already
#         carries the "{genome}" wildcard) -- not fixed_path()/db_token(),
#         which only know about the single-genome get_workflow_prefix().
# ────────────────────────────────────────────────────────────────────────────
# MYTOOL_OUTPUT = f"{GENOME_PREFIX}mytool/results.tsv"
# MYTOOL_TOKEN = f"{GENOME_PREFIX}mytool/mytool_db.tkn"
#
# STEP 2: Add to rule all, expanded over every discovered genome, gated on
#         an opt-out config flag (defaults True so nothing changes until a
#         caller actually sends run_mytool=false -- this is what would let a
#         future frontend checkbox disable a tool per run).
# ────────────────────────────────────────────────────────────────────────────
# rule all:
#     input:
#         ...,
#         (expand(MYTOOL_TOKEN, genome=list(GENOMES.keys())) if rc_bool('run_mytool', True, config=config) else [])
#
# STEP 3: Copy and customize this template. input: is a lambda so it can
# look up THIS genome's real file via wildcards.genome -- swap
# GENOMES[wildcards.genome] for whatever upstream rule's output this tool
# actually needs (e.g. another tool's per-genome .faa).
# ────────────────────────────────────────────────────────────────────────────
#
# rule run_MYTOOL:
#     """Brief description of what MYTOOL does"""
#     input:
#         lambda wildcards: GENOMES[wildcards.genome]
#     output:
#         results=MYTOOL_OUTPUT
#     group: "MYTOOL"
#     threads: rc('MYTOOL.threads', 4, config=config)
#     resources:
#         mem_mb=rc('MYTOOL.mem_mb', 4000, config=config),
#         runtime=rc('MYTOOL.runtime', 120, config=config)
#     params:
#         output_dir=lambda wildcards: rc('MYTOOL.output_dir', f'mytool_work/{wildcards.genome}', config=config),
#         db=db_path('MYTOOL', config=config, workflow_id='margie_sb')
#     container: sif_path('MYTOOL.sif', config=config, workflow_id='margie_sb')
#     shell:
#         """
#         echo "=== MARGIE_SB PHASE N: MYTOOL ({wildcards.genome}) ==="
#         /usr/local/bin/run -i {input} -o {params.output_dir} -d {params.db} -t {threads}
#         cp $(find {params.output_dir} -name "results.tsv") {output.results}
#         """
#
#
# rule load_MYTOOL_to_db:
#     """Load MYTOOL results into SQLite database"""
#     input:
#         results=MYTOOL_OUTPUT
#     output:
#         tkn=MYTOOL_TOKEN
#     group: "MYTOOL"
#     params:
#         db=MAIN_DATABASE,
#         script=LOAD_SCRIPT
#     shell:
#         """
#         {LOADER_PYTHON} {params.script} tsv {input.results} {params.db} mytool --token {output.tkn}
#         """
#
# ─────────────────────────────────────────────────────────────────────────────
# NOTES:
# ─────────────────────────────────────────────────────────────────────────────
# • Build every path from GENOME_PREFIX (carries the {genome} wildcard), not
#   fixed_path()/build_filepath()/db_token() -- those only know the single
#   global input_fasta, not a per-genome wildcard.
# • Use a per-genome params.output_dir (f'..._work/{wildcards.genome}') for
#   any tool's scratch/working directory -- never a bare relative constant,
#   or two genomes running in parallel will collide in the same folder.
# • Container path: sif_path('TOOL.sif', config=config, workflow_id='margie_sb').
#   DB path (if the tool needs one): db_path('tool_key', config=config, workflow_id='margie_sb').
# • Always use rc() for configurable parameters with sensible defaults.
# • Group name should match tool name for easier debugging.
# • loader format: tsv, csv, gff (see load_to_db.py for supported formats).
# • Start each run_MYTOOL rule's shell: with an
#   echo "=== MARGIE_SB PHASE N: MYTOOL ({wildcards.genome}) ===" line. "Phase"
#   here is just a conceptual ordering label, not tracked anywhere in code --
#   but the frontend's job page displays raw stdout/stderr verbatim in its
#   Logs panel, so this one line gives every run a readable, consistent,
#   phase-by-phase trail with zero backend/API/frontend changes needed.
# • If a tool's container needs no database (e.g. quast), drop params.db and
#   the -d flag -- not every MYTOOL invocation needs every template line.
# • If run_MYTOOL needs a non-default resources.slurm_partition (e.g.
#   gtdbtk needs 'highmem'), do NOT give its paired load_MYTOOL_to_db rule
#   the same group: -- the SLURM executor submits a shared group as ONE
#   job, and a single job can't request two different partitions. Only
#   share a group: between rules that need the exact same resources.
# ═════════════════════════════════════════════════════════════════════════════