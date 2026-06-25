"""
Shared utilities for Snakemake workflows.

Provides configuration helpers and standardized path generation.
Import into Snakemake files to access config and generate output paths.
"""
from pathlib import Path

# Per-workflow fallback defaults for the root-path settings (sif_path, db_root,
# input_path, output_path). Single source of truth — consulted by
# workflow_registry.workflow_path_params() (what the Profile UI shows/saves)
# and by sif_dir()/db_path()/do_margie_sb()'s own fallback chain, so a
# workflow's defaults only need to be written in one place.
WORKFLOW_PATH_DEFAULTS: dict[str, dict[str, str]] = {
    'margie_sb': {
        'sif_path': '/scratch/negishi/bhattar3/margie/sif',
        'db_root': '/scratch/negishi/bhattar3/margie/db',
        'input_path': '/scratch/negishi/bhattar3/margie/user-input',
        'output_path': '/scratch/negishi/bhattar3/margie/output',
    },
    'margie': {
        # No input_path default: margie takes one specific genome file, not
        # a folder, so there's no generally-correct file to default to.
        'output_path': '/scratch/negishi/bhattar3/margie/output',
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
    """Like rc(), but for boolean flags (e.g. 'run_quast') -- handles both a
    real Python bool (set via a YAML --configfile, which parses true/false
    natively) and a string (set via `--config key=value` on the CLI, which
    Snakemake always keeps as a string -- so a bare `if rc(...)` would treat
    the string 'false' as truthy and silently never disable anything)."""
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
    """Expand user-home markers for config-supplied filesystem paths."""
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
    """Same as get_workflow_prefix(), but for an explicit genome stem rather
    than the single global input_fasta. Used for multi-genome batch runs
    where each genome gets its own '{output_dir}/{genome_stem}/' subfolder.
    """
    output_dir = config.get('output_dir', '')
    if not output_dir:
        return ''

    abs_path = Path(output_dir).resolve()
    return f"{abs_path}/{genome_stem}/"


def get_container_outputs_prefix_for(genome_stem: str, config) -> str:
    """Per-genome root for each tool's full container output (raw/processed/
    pipeline-log/etc), separate from get_workflow_prefix_for()'s output_dir,
    which Snakemake tracks as real rule outputs and stays lean. Always
    '{output_dir}/original_container_outputs/{genome_stem}/'."""
    output_dir = config.get('output_dir', '')
    if not output_dir:
        return ''

    abs_path = Path(output_dir).resolve() / 'original_container_outputs'
    return f"{abs_path}/{genome_stem}/"


GENOME_EXTENSIONS = ('.fasta', '.fa', '.fna', '.fasta.gz', '.fa.gz', '.fna.gz')


def discover_genomes(input_path: str, recursive: bool = False) -> dict[str, str]:
    """Resolve an input path to a {stem: filepath} map. A single file
    resolves to one entry; a directory is scanned for recognized genome
    files, one entry per file, sorted by name.

    Non-recursive by default, so nested reference-genome folders (e.g.
    margie_sb's synteny-input/<genome>/...) aren't swept into the primary
    genome set. recursive=True opts every nested genome in as a primary
    target instead of a reference."""
    path = Path(config_path(input_path))
    if path.is_dir():
        entries = path.rglob('*') if recursive else path.iterdir()
        genomes = {}
        for entry in sorted(entries):
            if entry.is_file() and entry.name.lower().endswith(GENOME_EXTENSIONS):
                genomes[entry.stem] = str(entry)
        return genomes
    return {path.stem: str(path)}


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
    """Return the configured SIF directory, expanded to an absolute user path.

    Looks up '<workflow_id>.sif_path' when workflow_id is given (the
    per-workflow setting set in the Profile page's "Workflow Specific
    Settings"), otherwise falls back to a bare 'sif_path' key.
    """
    if config is None:
        raise ValueError("config parameter is required")

    key = f'{workflow_id}.sif_path' if workflow_id else 'sif_path'
    fallback = WORKFLOW_PATH_DEFAULTS.get(workflow_id, {}).get('sif_path', '~/.cache/bioinformatics-tools')
    return config_path(rc(key, fallback, config=config))


def sif_path(filename: str, config=None, workflow_id: str | None = None) -> str:
    """Resolve a container filename relative to the configured SIF directory."""
    if config is None:
        raise ValueError("config parameter is required")

    return str(Path(sif_dir(config=config, workflow_id=workflow_id)) / filename)


def db_path(tool: str, config=None, default_root: str | None = None,
            workflow_id: str | None = None) -> str:
    """Resolve a tool database path from either db.<tool> or <workflow_id>.db_root/<tool>.

    Not every tool needs a database — this just resolves where to *look* for
    one when a tool does need it.
    """
    if config is None:
        raise ValueError("config parameter is required")

    explicit = rc(f'db.{tool}', None, config=config)
    if explicit:
        return config_path(explicit)

    if default_root is None:
        default_root = WORKFLOW_PATH_DEFAULTS.get(workflow_id, {}).get('db_root', '/depot/lindems/data/Databases')

    db_root_key = f'{workflow_id}.db_root' if workflow_id else 'db_root'
    db_root = rc(db_root_key, default_root, config=config)
    return str(Path(config_path(db_root)) / tool)
