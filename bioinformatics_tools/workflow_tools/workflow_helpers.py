"""
Shared utilities for Snakemake workflows.

Provides configuration helpers and standardized path generation.
Import into Snakemake files to access config and generate output paths.
"""
import logging
from pathlib import Path

LOGGER = logging.getLogger(__name__)

# Repo root (the bioinformatics-tools checkout), used for machine-agnostic fallback paths.
_REPO_ROOT = Path(__file__).resolve().parents[2]

# Per-workflow fallback defaults for the root-path settings, read by
# workflow_registry.workflow_path_params() (the Profile UI) and by
# sif_dir()/db_path()/do_margie_sb(). Input/output paths stay user-provided.
WORKFLOW_PATH_DEFAULTS: dict[str, dict[str, str]] = {
    'margie_sb': {
        'sif_path': '',
        # The tools' reference databases (read only, shared by everyone).
        'db_root': '',
        'input_path': '',
        'output_path': '',
    },
    'margie': {
        # No input_path default: margie takes one specific genome file.
        'output_path': '',
    },
}


def rc(key: str, default:str|None = None, config=None):
    """
    Rule Config: Get config value using dot notation for nested keys.

    Supports arbitrary nesting via dot notation:
    key: Config key using dot notation (e.g., 'prodigal.mem_mb')

    Returns:
        The config value OR default OR if not set
    """
    if config is None:
        raise ValueError("config parameter is required")

    # Split key by dots and traverse nested dict
    parts = key.split('.')
    value = config

    for part in parts:
        if isinstance(value, dict) and part in value:
            value = value[part]
        else:
            return default

    return value


def rc_bool(key: str, default: bool = True, config=None) -> bool:
    """Like rc(), for boolean flags such as 'run_quast'.

    Accepts a real bool (YAML --configfile) or a string (`--config key=value`),
    so the string 'false' is not treated as truthy."""
    value = rc(key, default, config=config)
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        return value.strip().lower() not in ('false', '0', 'no', 'off', '')
    return bool(value)


def get_input_stem(config):
    """
    Extract stem (filename without extension) from input file.
    """
    input_file = config.get('input_fasta', 'output')
    return Path(input_file).stem


def config_path(value: str | Path | None) -> str:
    """Expands user-home markers in config-supplied filesystem paths."""
    if value is None:
        return ''
    return str(Path(value).expanduser())


def get_workflow_prefix(config) -> str | Path:
    """
    Get output directory prefix with stem subdirectory and trailing slash.

    Each genome gets isolated in its own subdirectory for clean batch processing.

    Example:
        input_fasta: 'ecoli.fasta'
        output_dir: '/results/batch-run'
        → '/results/batch-run/ecoli/'
    """
    output_dir = config.get('output_dir', '')
    if not output_dir:
        return ''

    stem = get_input_stem(config)
    abs_path = Path(output_dir).resolve()
    return f"{abs_path}/{stem}/"


def get_workflow_prefix_for(genome_stem: str, config) -> str:
    """Returns get_workflow_prefix() for an explicit genome stem, i.e. that
    genome's '{output_dir}/{genome_stem}/' folder in multi-genome runs.
    """
    output_dir = config.get('output_dir', '')
    if not output_dir:
        return ''

    abs_path = Path(output_dir).resolve()
    return f"{abs_path}/{genome_stem}/"


def get_container_outputs_prefix_for(genome_stem: str, config) -> str:
    """Returns the per-genome root for each tool's full container output,
    '{output_dir}/original_container_outputs/{genome_stem}/', kept apart from
    the lean rule outputs Snakemake tracks."""
    output_dir = config.get('output_dir', '')
    if not output_dir:
        return ''

    abs_path = Path(output_dir).resolve() / 'original_container_outputs'
    return f"{abs_path}/{genome_stem}/"


GENOME_EXTENSIONS = ('.fasta', '.fa', '.fna', '.fasta.gz', '.fa.gz', '.fna.gz')


def genome_stem(name: str) -> str:
    """Returns a genome file's key: the basename minus .gz and its recognised extension.

    .fasta and .fna map to the same key.
    """
    lowered = name.lower()
    if lowered.endswith('.gz'):
        name, lowered = name[:-3], lowered[:-3]
    for ext in ('.fasta', '.fna', '.fa'):
        if lowered.endswith(ext):
            return name[: -len(ext)]
    return Path(name).stem


def default_store_root() -> str:
    """Returns the default store root, /scratch/<cluster>/<user>/margie-2026.

    A fallback for runs started outside the app, which always names the stores
    (api/services/user_stores.py); depot base copies are never written."""
    import os
    import socket
    cluster = (socket.getfqdn().split('.') + ['', ''])[1] or 'cluster'
    return f"/scratch/{cluster}/{os.environ.get('USER', 'user')}/margie-2026"


def genome_calls(genomes: dict[str, str], config) -> dict[str, dict[str, str]]:
    """Returns each genome's domain, genetic code and gene caller from config alone.

    `margie_sb.genome_info` in config.yaml maps a genome file name or stem to
    {domain, genetic_code}. RASTtk needs both; a genome missing either is called
    by Prodigal unless GTDB-Tk runs, which determines both. Computed at parse
    time so the Snakefile can split genomes between run_rasttk and run_prodigal.
    """
    run_gtdbtk = rc_bool('run_gtdbtk', True, config=config)
    meta = rc('margie_sb.genome_info', None, config=config) or rc('genome_info', None, config=config) or {}
    if not isinstance(meta, dict):
        meta = {}
    by_key = {}
    for key, row in meta.items():
        if isinstance(row, dict):
            by_key.setdefault(str(key), row)
            by_key.setdefault(genome_stem(str(key)), row)
    out = {}
    for genome, fasta in genomes.items():
        row = by_key.get(Path(fasta).name) or by_key.get(genome) or {}
        raw = str(row.get('domain') or '').strip().lower()
        domain = 'Bacteria' if raw.startswith('b') else 'Archaea' if raw.startswith('a') else ''
        code = str(row.get('genetic_code') or '').strip()
        code = code if code.isdigit() else ''
        caller = 'rasttk' if run_gtdbtk or (domain and code) else 'prodigal'
        out[genome] = {'domain': domain, 'genetic_code': code, 'gene_caller': caller,
                       'source': 'metadata' if domain or code else ('gtdbtk' if run_gtdbtk else 'none')}
    return out


def discover_genomes(input_path: str, recursive: bool = False) -> dict[str, str]:
    """Resolves an input path to a {stem: filepath} map.

    A file gives one entry; a directory is scanned for recognised genome files,
    sorted by name. Non-recursive by default, so nested reference genomes (e.g.
    synteny-input/<genome>/...) stay out; recursive=True includes them."""
    path = Path(config_path(input_path))
    if path.is_dir():
        entries = path.rglob('*') if recursive else path.iterdir()
        genomes = {}
        for entry in sorted(entries):
            if entry.is_file() and entry.name.lower().endswith(GENOME_EXTENSIONS):
                stem = genome_stem(entry.name)
                previous = genomes.get(stem)
                if previous is not None:
                    # Two files sharing a stem are the same genome if their content
                    # hash matches (either is kept); otherwise it is an error. Hashing
                    # happens only on a stem collision.
                    from bioinformatics_tools.workflow_tools.load_to_db import compute_fasta_hash
                    if compute_fasta_hash(previous) == compute_fasta_hash(str(entry)):
                        LOGGER.info(
                            "Genome '%s' present as both %s and %s with identical "
                            "content -- same genome, using %s",
                            stem, Path(previous).name, entry.name, Path(previous).name)
                        continue
                    raise ValueError(
                        f"Two DIFFERENT genomes map to the same name '{stem}': "
                        f"{previous} and {entry}. Their contents differ, so one "
                        f"would overwrite the other's outputs. Rename one so each "
                        f"genome has a distinct name.")
                genomes[stem] = str(entry)
        return genomes
    return {genome_stem(path.name): str(path)}


def fixed_path(relative_path: str, config=None) -> str:
    """Prepend workflow prefix to a relative path. E.g., 'prodigal/file.faa' → 'results/prodigal/file.faa'"""
    if config is None:
        raise ValueError("config parameter is required")

    prefix = get_workflow_prefix(config)
    return f"{prefix}{relative_path}"


def build_filepath(config_string: str, suffix: str, default: str = None, config=None) -> str:
    """
    Build filepath with config lookup and auto-generation.

    Priority: 1) Config value, 2) Default path, 3) Auto-generate {tool}/{stem}-{tool}.{suffix}

    Examples:
        build_filepath('pfam.output', suffix='tsv') → 'results/pfam/ecoli-pfam.tsv'
        build_filepath('cog.output', suffix='txt', default='cog/results.txt') → 'results/cog/results.txt'
        # With config pfam.output='custom.tsv' → 'results/pfam/custom.tsv'
    """
    if config is None:
        raise ValueError("config parameter is required")

    prefix = get_workflow_prefix(config)

    parts = config_string.split('.')
    if len(parts) < 1:
        raise ValueError(f"config_string must have at least one part, got: {config_string}")
    tool = parts[0]

    # Check config first
    config_filename = rc(config_string, None, config=config)
    if config_filename:
        return f"{prefix}{tool}/{config_filename}"

    # Use default if provided
    if default:
        return f"{prefix}{default}"

    # Auto-generate with stem
    stem = get_input_stem(config)
    return f"{prefix}{tool}/{stem}-{tool}.{suffix}"


def db_token(tool, config=None):
    """Generate database token path: {tool}/{tool}_db.tkn"""
    return fixed_path(f'{tool}/{tool}_db.tkn', config=config)


def sif_dir(config=None, workflow_id: str | None = None) -> str:
    """Returns the configured SIF directory as an absolute user path.

    Reads '<workflow_id>.sif_path' when workflow_id is given, else 'sif_path'.
    """
    if config is None:
        raise ValueError("config parameter is required")

    key = f'{workflow_id}.sif_path' if workflow_id else 'sif_path'
    fallback = WORKFLOW_PATH_DEFAULTS.get(workflow_id, {}).get('sif_path', '~/.cache/bioinformatics-tools')
    return config_path(rc(key, fallback, config=config))


def sif_path(filename: str, config=None, workflow_id: str | None = None) -> str:
    """Resolves a container filename relative to the configured SIF directory."""
    if config is None:
        raise ValueError("config parameter is required")

    return str(Path(sif_dir(config=config, workflow_id=workflow_id)) / filename)


def db_path(tool: str, config=None, default_root: str | None = None,
            workflow_id: str | None = None) -> str:
    """Resolves a tool database path from db.<tool> or <workflow_id>.db_root/<tool>."""
    if config is None:
        raise ValueError("config parameter is required")

    explicit = rc(f'db.{tool}', None, config=config)
    if explicit:
        return config_path(explicit)

    if default_root is None:
        default_root = WORKFLOW_PATH_DEFAULTS.get(workflow_id, {}).get('db_root', '')

    db_root_key = f'{workflow_id}.db_root' if workflow_id else 'db_root'
    db_root = rc(db_root_key, default_root, config=config)
    return str(Path(config_path(db_root)) / tool)
