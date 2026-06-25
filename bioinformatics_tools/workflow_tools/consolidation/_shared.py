#!/usr/bin/env python3
"""Shared discovery/normalisation/row-key logic for margie_sb's
consolidation pipeline (detect-columns.py, merge-all-columns.py, and
later filter-for-labeling.py all import from here, instead of each
duplicating the same ~150 lines).

Not meant to be run directly -- leading underscore signals "internal".
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

# Identity-like columns pulled directly from rasttk's raw per-row values
# in build_merged_rows() (not via merge_rows_for_gene(), which only knows
# about tool_columns -- these are GLOBAL_COLUMNS, deliberately excluded
# from that set so they don't get double-prefixed or multi-hit-wrapped).
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

# Excluded from the merge entirely: meta-stages, comparative/genome-level
# tools, genome-level QC/taxonomy tools with no feature_id, and envelope
# (its decision is already propagated into deepsig/psortb/signalp4's own
# results.tsv via enrich_with_envelope.py -- see GLOBAL_COLUMNS above).
ALWAYS_EXCLUDED: frozenset[str] = frozenset({
    "consolidation", "labeling", "fingerprint", "scoring", "scoring_heuristic",
    "fingerprint_database", "synteny", "aai", "ani", "closest_organisms",
    "gtdbtk", "quast", "envelope",
})

SPECIAL_FILENAME_OVERRIDES: dict[str, str] = {
    "rasttk": "rast.tsv",
}


# ─── Per-tool row-key strategy ─────────────────────────────────────────────────
#
# Determines how multiple rows for the SAME gene within ONE tool's table get
# joined into a single cell. Only used when a gene actually has >1 row for
# that tool -- a single-row gene's values stay bare, no key: prefix at all.

class KeySpec(NamedTuple):
    kind: str                    # "column" | "range" | "composite" | "none" | "synthetic_index"
    columns: tuple[str, ...]     # column name(s) feeding the key, post-normalise_col


def _accession_key(prefix: str) -> KeySpec:
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
    # geneprop: keyed by GENEPROP_id itself (its real per-row grain is
    # finer -- one row per (tigrfam_hit, genprop_id, step) -- but the step
    # level isn't useful in the merged table; collapsing to one entry per
    # distinct GENEPROP_id, with GENEPROP_description specially rendered
    # as "id: description(status)" -- see merge-all-columns.py).
    "geneprop": _accession_key("GENEPROP"),
    # Always exactly one row per gene -- no wrapping ever needed.
    "psortb": KeySpec("none", ()), "signalp4": KeySpec("none", ()),
    "signalp6": KeySpec("none", ()), "operon": KeySpec("none", ()),
    "rasttk": KeySpec("none", ()),
}

# Columns that should always render as their own bare deduplicated value
# (like the key column itself), never wrapped in "key: value" -- e.g.
# geneprop's tigrfam_hit is effectively a second, constant-per-group key
# ("which TIGRFAM triggered this"), not a per-GenProp varying value.
BARE_VALUE_COLUMNS: dict[str, set[str]] = {
    "geneprop": {"GENEPROP_tigrfam_hit"},
}


def resolve_key_spec(tool_name: str) -> KeySpec:
    if tool_name in ROW_KEY_SPEC:
        return ROW_KEY_SPEC[tool_name]
    if tool_name.startswith("interpro_"):
        db = tool_name[len("interpro_"):]
        return _accession_key(f"INTERPRO_{db.upper()}")
    # Unknown future tool: fall back to a synthetic per-row index rather
    # than crash -- still correct, just less informative than a real key.
    return KeySpec("synthetic_index", ())


def format_composite_key(values: dict[str, str], spec: KeySpec) -> str:
    if spec.kind == "column":
        return values.get(spec.columns[0], "")
    if spec.kind == "range":
        start, end = (values.get(c, "") for c in spec.columns)
        return f"{start}-{end}"
    return ""


# ─── Helpers ──────────────────────────────────────────────────────────────────

def feature_id_invalid(fid: str) -> bool:
    if not fid or len(fid) < 3 or " " in fid:
        return True
    if fid.startswith(("5'", "3'", "5`", "3`")):
        return True
    return False


def normalise_col(raw_col: str, tool_name: str) -> str:
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
    """Find every tool's results.tsv directly under output_root/<tool>/.

    margie_sb's output/ tree is flat -- output_root/<tool>/<tool>_results.tsv,
    no per-organism/processed/ nesting (that nesting only exists in the
    separate container_outputs/ tree, which holds raw working data, not
    final results). interpro is the one exception with multiple files, one
    per member database actually run; rasttk's real filename is rast.tsv.
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
    tool_name: str
    source_path: Path
    rows_by_feature: dict[str, list[dict[str, str]]]
    tool_columns: list[str]          # canonicalized (post normalise_col), in load order
    raw_columns: list[str]           # exactly as they appeared in the source file


def load_tool_table(tool_name: str, source_path: Path) -> ToolTable:
    rows_by_feature: dict[str, list[dict[str, str]]] = defaultdict(list)
    tool_columns: list[str] = []
    seen_cols: set[str] = set()
    raw_columns: list[str] = []

    with open(source_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        raw_columns = list(reader.fieldnames or [])

        # Register columns from the HEADER itself, not just rows that
        # successfully load -- a tool with zero data rows for this genome
        # (e.g. dbcan finding nothing) would otherwise contribute zero
        # tool_columns, silently dropping DBCAN_id/DBCAN_description from
        # every downstream merged/filtered table instead of leaving them
        # present-but-empty.
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
#
# Shared by every script that reads merge-all-columns.py's output and needs
# to know "did this tool find anything for this gene" -- extract-hit-counts.py
# and extract-database-coverage.py both import these instead of duplicating.

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
    if not value:
        return 0
    return len([v for v in value.split(";") if v])


def count_semicolon_space_list(value: str) -> int:
    if not value:
        return 0
    return len([v for v in value.split("; ") if v])


def tool_hit_count(row: dict[str, str], tool: str) -> int:
    """Number of distinct hits `tool` has for this row, using whichever
    counting strategy matches its shape (see ID_BASED_TOOLS/
    SEGMENT_BASED_TOOLS/SINGLE_ROW_TOOLS above)."""
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
