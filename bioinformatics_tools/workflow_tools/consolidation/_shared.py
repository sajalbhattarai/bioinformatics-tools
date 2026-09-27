#!/usr/bin/env python3
"""Shared tool-table discovery, column normalisation and row-key logic for
the consolidation scripts (detect-columns.py, merge-all-columns.py and the
derived-table scripts import it). Not run directly.
"""
from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path
from typing import NamedTuple

csv.field_size_limit(10_000_000)

# ─── Column aliases ────────────────────────────────────────────────────────────
COLUMN_ALIAS_MAP: dict[str, str] = {
    "organism": "organism_name",   # RAST
}

COLUMNS_DROPPED: frozenset[str] = frozenset()

GLOBAL_COLUMNS: frozenset[str] = frozenset({
    "organism_name", "domain", "feature_id",
    "gene_id", "gene_start", "gene_end", "na_length", "aa_length",
    "na_seq", "aa_seq",
    "ENVELOPE_envelope_type", "ENVELOPE_inference_basis", "ENVELOPE_evidence_json",
})

# Identity columns copied straight from rasttk's rows in build_merged_rows();
# they are GLOBAL_COLUMNS, so they are never prefixed or multi-hit-wrapped.
RASTTK_IDENTITY_COLUMNS: tuple[str, ...] = (
    "gene_id", "gene_start", "gene_end", "na_length", "aa_length", "na_seq", "aa_seq",
)

TOOL_PREFIX_MAP: dict[str, str] = {
    "pfam_": "PFAM_", "tigrfam_": "TIGRFAM_", "pgap_": "PGAP_",
    "deepsig_": "DEEPSIG_", "tmbed_": "TMBED_", "psortb_": "PSORTB_",
    "geneprop_": "GENEPROP_", "dbcan_": "DBCAN_", "kegg_": "KEGG_",
    "cog_": "COG_", "merops_": "MEROPS_", "tcdb_": "TCDB_",
    "uniprot_": "UNIPROT_", "eggnog_": "EGGNOG_", "interpro_": "INTERPRO_",
    "rast_": "RAST_", "rasttk_": "RASTTK_", "signalp4_": "SIGNALP4_",
    "signalp6_": "SIGNALP6_", "phobius_": "PHOBIUS_", "operon_": "OPERON_",
}

# Stages excluded from the merge: meta-stages, genome-level tools without a
# feature_id, and envelope (already folded into deepsig/psortb/signalp4 by
# enrich_with_envelope.py).
ALWAYS_EXCLUDED: frozenset[str] = frozenset({
    "consolidation", "labeling", "fingerprint", "scoring", "scoring_heuristic",
    "fingerprint_database", "synteny", "aai", "ani", "closest_organisms",
    "gtdbtk", "quast", "envelope",
})

SPECIAL_FILENAME_OVERRIDES: dict[str, str] = {
    "rasttk": "rast.tsv",
}


# ─── Per-tool row-key strategy ─────────────────────────────────────────────────
# Decides how several rows for one gene in one tool's table join into a cell;
# a gene with a single row keeps bare values without a "key:" prefix.

class KeySpec(NamedTuple):
    """Describes how a tool's per-gene row key is formed."""
    kind: str                    # "column" | "range" | "composite" | "none" | "synthetic_index"
    columns: tuple[str, ...]     # column name(s) feeding the key, post-normalise_col


def _accession_key(prefix: str) -> KeySpec:
    """Returns a KeySpec keyed on the tool's <PREFIX>_id column."""
    return KeySpec("column", (f"{prefix}_id",))


ROW_KEY_SPEC: dict[str, KeySpec] = {
    # Accession/family-ID tools: the id column itself is the key.
    "cog": _accession_key("COG"), "eggnog": _accession_key("EGGNOG"),
    "kegg": _accession_key("KEGG"), "pfam": _accession_key("PFAM"),
    "tigrfam": _accession_key("TIGRFAM"), "pgap": _accession_key("PGAP"),
    "merops": _accession_key("MEROPS"), "tcdb": _accession_key("TCDB"),
    "dbcan": _accession_key("DBCAN"), "uniprot": _accession_key("UNIPROT"),
    # Coordinate-range tools: no accession exists, but start/end does.
    "phobius": KeySpec("range", ("PHOBIUS_segment_start", "PHOBIUS_segment_end")),
    "tmbed":   KeySpec("range", ("TMBED_segment_start", "TMBED_segment_end")),
    "deepsig": KeySpec("range", ("DEEPSIG_start", "DEEPSIG_end")),
    # geneprop: collapses its (tigrfam_hit, genprop_id, step) rows to one entry
    # per GENEPROP_id; merge-all-columns.py renders the description specially.
    "geneprop": _accession_key("GENEPROP"),
    # Tools with exactly one row per gene.
    "psortb": KeySpec("none", ()), "signalp4": KeySpec("none", ()),
    "signalp6": KeySpec("none", ()), "operon": KeySpec("none", ()),
    "rasttk": KeySpec("none", ()),
}

# Columns rendered as bare deduplicated values, never "key: value"
# (geneprop's tigrfam_hit is constant per group, like a second key).
BARE_VALUE_COLUMNS: dict[str, set[str]] = {
    "geneprop": {"GENEPROP_tigrfam_hit"},
}


def resolve_key_spec(tool_name: str) -> KeySpec:
    """Returns the row-key strategy for a tool, deriving InterPro member keys by name."""
    if tool_name in ROW_KEY_SPEC:
        return ROW_KEY_SPEC[tool_name]
    if tool_name.startswith("interpro_"):
        db = tool_name[len("interpro_"):]
        return _accession_key(f"INTERPRO_{db.upper()}")
    # Unknown tools fall back to a synthetic per-row index.
    return KeySpec("synthetic_index", ())


def format_composite_key(values: dict[str, str], spec: KeySpec) -> str:
    """Builds the key string for one row: the id column value or "start-end"."""
    if spec.kind == "column":
        return values.get(spec.columns[0], "")
    if spec.kind == "range":
        start, end = (values.get(c, "") for c in spec.columns)
        return f"{start}-{end}"
    return ""


# ─── Helpers ──────────────────────────────────────────────────────────────────

def feature_id_invalid(fid: str) -> bool:
    """Returns True for empty, short, spaced or 5'/3'-prefixed feature ids."""
    if not fid or len(fid) < 3 or " " in fid:
        return True
    if fid.startswith(("5'", "3'", "5`", "3`")):
        return True
    return False


def normalise_col(raw_col: str, tool_name: str) -> str:
    """Maps a raw column name to its canonical TOOL_-prefixed form via TOOL_PREFIX_MAP."""
    col = COLUMN_ALIAS_MAP.get(raw_col, raw_col)
    if col in GLOBAL_COLUMNS:
        return col
    col_lower = col.lower()
    for prefix_lower, prefix_upper in TOOL_PREFIX_MAP.items():
        if col_lower.startswith(prefix_lower):
            return prefix_upper + col[len(prefix_lower):]
    tool_lower = tool_name.lower()
    if col_lower.startswith(tool_lower + "_"):
        return tool_name.upper() + col[len(tool_name):]
    return f"{tool_name.upper()}_{col}"


# ─── Discovery ────────────────────────────────────────────────────────────────

def discover_tool_tables(output_root: Path, extra_excluded: set[str]) -> list[tuple[str, Path]]:
    """Lists (tool_name, path) for each tool's results table under output_root/<tool>/.

    InterPro yields one entry per member-database file; rasttk's file is rast.tsv.
    """
    excluded = ALWAYS_EXCLUDED | extra_excluded
    pairs: list[tuple[str, Path]] = []
    if not output_root.is_dir():
        return pairs

    for tool_dir in sorted(output_root.iterdir()):
        if not tool_dir.is_dir():
            continue
        tool_name = tool_dir.name
        if tool_name in excluded:
            continue

        if tool_name == "interpro":
            for tsv in sorted(tool_dir.glob("interpro_*_results.tsv")):
                if tsv.name == "interpro_results.tsv":
                    continue  # skip the unified (row-per-hit) file
                sub = tsv.name[len("interpro_"): -len("_results.tsv")]
                pairs.append((f"interpro_{sub}", tsv))
            continue

        filename = SPECIAL_FILENAME_OVERRIDES.get(tool_name, f"{tool_name}_results.tsv")
        canonical = tool_dir / filename
        if canonical.is_file():
            pairs.append((tool_name, canonical))
            continue

        for tsv in sorted(tool_dir.glob("*.tsv")):
            pairs.append((tool_name, tsv))
            break

    return pairs


# ─── Loading ──────────────────────────────────────────────────────────────────

class ToolTable(NamedTuple):
    """One loaded tool table: rows grouped by feature_id plus its column lists."""
    tool_name: str
    source_path: Path
    rows_by_feature: dict[str, list[dict[str, str]]]
    tool_columns: list[str]          # canonicalized (post normalise_col), in load order
    raw_columns: list[str]           # exactly as they appeared in the source file


def load_tool_table(tool_name: str, source_path: Path) -> ToolTable:
    """Reads a tool's TSV with csv.DictReader, normalising columns and grouping rows by feature_id."""
    rows_by_feature: dict[str, list[dict[str, str]]] = defaultdict(list)
    tool_columns: list[str] = []
    seen_cols: set[str] = set()
    raw_columns: list[str] = []

    with open(source_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        raw_columns = list(reader.fieldnames or [])

        # Registers columns from the header so a tool with no rows still
        # contributes its (empty) columns downstream.
        for raw_col in raw_columns:
            if raw_col in COLUMNS_DROPPED:
                continue
            ncol = normalise_col(raw_col, tool_name)
            if ncol not in GLOBAL_COLUMNS and ncol not in seen_cols:
                seen_cols.add(ncol)
                tool_columns.append(ncol)

        for raw_row in reader:
            norm: dict[str, str] = {}
            fid = ""
            for raw_col, val in raw_row.items():
                val = (val or "").strip()
                if raw_col in COLUMNS_DROPPED:
                    continue
                ncol = normalise_col(raw_col, tool_name)
                if ncol == "feature_id":
                    fid = val
                norm[ncol] = val
                if ncol not in GLOBAL_COLUMNS and ncol not in seen_cols:
                    seen_cols.add(ncol)
                    tool_columns.append(ncol)
            if feature_id_invalid(fid):
                continue
            rows_by_feature[fid].append(norm)

    return ToolTable(tool_name, source_path, dict(rows_by_feature), tool_columns, raw_columns)


# ─── Derived-table tool categorization ─────────────────────────────────────────
# Used by scripts that read merge-all-columns.py's output to count per-tool hits.

IDENTITY_COLUMNS: list[str] = [
    "feature_id", "organism_name", "domain",
    "gene_id", "gene_start", "gene_end",
    "na_length", "aa_length",
    "RAST_feature_type", "RAST_strand", "RAST_description",
    "na_seq", "aa_seq",
]

# tool_name -> source column carrying its bare ";"-joined accession list.
ID_BASED_TOOLS: dict[str, str] = {
    "COG": "COG_id", "KEGG": "KEGG_id", "EGGNOG": "EGGNOG_id",
    "PFAM": "PFAM_id", "TIGRFAM": "TIGRFAM_id", "PGAP": "PGAP_id",
    "MEROPS": "MEROPS_id", "TCDB": "TCDB_id", "DBCAN": "DBCAN_id",
    "UNIPROT": "UNIPROT_id", "GENEPROP": "GENEPROP_id",
    "INTERPRO": "INTERPRO_id",
    "INTERPRO_HAMAP": "INTERPRO_HAMAP_id", "INTERPRO_NCBIFAM": "INTERPRO_NCBIFAM_id",
    "INTERPRO_PFAM": "INTERPRO_PFAM_id", "INTERPRO_PANTHER": "INTERPRO_PANTHER_id",
    "INTERPRO_PIRSF": "INTERPRO_PIRSF_id", "INTERPRO_PIRSR": "INTERPRO_PIRSR_id",
    "INTERPRO_GENE3D": "INTERPRO_GENE3D_id", "INTERPRO_CDD": "INTERPRO_CDD_id",
    "INTERPRO_COILS": "INTERPRO_COILS_id", "INTERPRO_PRINTS": "INTERPRO_PRINTS_id",
    "INTERPRO_SMART": "INTERPRO_SMART_id", "INTERPRO_SFLD": "INTERPRO_SFLD_id",
    "INTERPRO_SUPERFAMILY": "INTERPRO_SUPERFAMILY_id", "INTERPRO_MOBIDB": "INTERPRO_MOBIDB_id",
    "INTERPRO_PROSITE_PATTERNS": "INTERPRO_PROSITE_PATTERNS_id",
    "INTERPRO_PROSITE_PROFILES": "INTERPRO_PROSITE_PROFILES_id",
    "INTERPRO_FUNFAM": "INTERPRO_FUNFAM_id", "INTERPRO_ANTIFAM": "INTERPRO_ANTIFAM_id",
}

# tool_name -> representative raw column, "; "-joined when multi-row.
SEGMENT_BASED_TOOLS: dict[str, str] = {
    "PHOBIUS": "PHOBIUS_segment_type",
    "TMBED": "TMBED_topology",
    "DEEPSIG": "DEEPSIG_feature_type",
}

# tool_name -> (column, value(s) counted as a positive call; None means
# "any non-empty value counts").
SINGLE_ROW_TOOLS: dict[str, tuple[str, frozenset[str] | None]] = {
    "PSORTB": ("PSORTB_localization", None),
    "SIGNALP4": ("SIGNALP4_is_signal_peptide", frozenset({"Y"})),
    "SIGNALP6": ("SIGNALP6_prediction", frozenset({"SP", "LIPO", "TAT", "TATLIPO", "PILIN"})),
    "OPERON": ("OPERON_id", None),
}


def count_id_list(value: str) -> int:
    """Counts non-empty entries in a ";"-joined list."""
    if not value:
        return 0
    return len([v for v in value.split(";") if v])


def count_semicolon_space_list(value: str) -> int:
    """Counts non-empty entries in a "; "-joined list."""
    if not value:
        return 0
    return len([v for v in value.split("; ") if v])


def tool_hit_count(row: dict[str, str], tool: str) -> int:
    """Returns the number of hits `tool` has in a merged row, by the tool's table shape."""
    if tool in ID_BASED_TOOLS:
        return count_id_list(row.get(ID_BASED_TOOLS[tool], ""))
    if tool in SEGMENT_BASED_TOOLS:
        return count_semicolon_space_list(row.get(SEGMENT_BASED_TOOLS[tool], ""))
    if tool in SINGLE_ROW_TOOLS:
        src_col, positive_values = SINGLE_ROW_TOOLS[tool]
        val = row.get(src_col, "")
        if positive_values is None:
            return 1 if val else 0
        return 1 if val in positive_values else 0
    raise KeyError(f"unknown tool: {tool}")


ALL_HIT_COUNT_TOOLS: list[str] = list(ID_BASED_TOOLS) + list(SEGMENT_BASED_TOOLS) + list(SINGLE_ROW_TOOLS)
