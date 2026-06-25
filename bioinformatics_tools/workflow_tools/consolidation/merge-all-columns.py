#!/usr/bin/env python3
"""merge-all-columns.py — Stage 2 of margie_sb's consolidation pipeline.

Merges every per-tool results.tsv under output/<genome>/ (rasttk through
every phase4-8 tool) into ONE wide per-feature (one row per gene) TSV with
every column from every tool, preserved exactly as it appears in each
tool's own table. This step ONLY merges -- it does not summarise, dedupe,
or pick "best" hits. The one thing it does do: when the *same gene* has
*multiple rows* within *one tool's* table (e.g. two PFAM domain hits, or
43 TMBED topology segments), those rows get joined into a single cell
using a "key: value; key2: value2" pairing, where "key" is whatever
already-existing column uniquely identifies each row for that tool (an
accession ID for database-hit tools, a start-end coordinate range for
segment/region tools, or a composite identifier for tools like geneprop
where the obvious ID column alone can repeat within one gene). Genes with
only ONE row for a given tool keep that value bare, with no key: prefix.

Discovery/normalisation/row-key logic lives in _shared.py (also used by
detect-columns.py). Downstream stages (filter-no-stat.py,
filter-for-labeling.py) read THIS script's output, not the original
per-tool files -- one expensive join, multiple cheap derived views.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from _shared import (
    BARE_VALUE_COLUMNS,
    RASTTK_IDENTITY_COLUMNS,
    ToolTable,
    discover_tool_tables,
    format_composite_key,
    load_tool_table,
    resolve_key_spec,
)

csv.field_size_limit(10_000_000)

# kegg_results.tsv stores EVERY KofamScan profile-vs-query comparison,
# significant or not, so a single gene can carry hundreds of rows where
# only a handful clear KEGG_is_above_threshold=true. Every other tool's
# own search command already filters at the source (e.g.
# pfam/tigrfam/pgap's --cut_ga/--cut_tc), so this filter is kegg-specific,
# not a general rule. Applied here (not in _shared.load_tool_table) so
# detect-columns.py still surfaces the raw, unfiltered reality.
TOOL_ROW_FILTERS: dict[str, tuple[str, str]] = {
    "kegg": ("KEGG_is_above_threshold", "true"),
}


def apply_tool_row_filter(table: ToolTable) -> ToolTable:
    filt = TOOL_ROW_FILTERS.get(table.tool_name)
    if not filt:
        return table
    col, want = filt
    filtered = {
        fid: kept
        for fid, rows in table.rows_by_feature.items()
        if (kept := [r for r in rows if r.get(col, "").lower() == want])
    }
    return table._replace(rows_by_feature=filtered)


# ─── Merge ────────────────────────────────────────────────────────────────────

def merge_rows_for_gene(tool: ToolTable, fid: str) -> dict[str, str]:
    """Merge one tool's row(s) for one gene into a single dict of columns.

    A single row's values pass through bare. Multiple rows get joined per
    column as "key1: val1; key2: val2", deduplicated by key (first
    occurrence wins), using that tool's row-key strategy.
    """
    rows = tool.rows_by_feature.get(fid, [])
    if not rows:
        return {col: "" for col in tool.tool_columns}
    if len(rows) == 1:
        return {col: rows[0].get(col, "") for col in tool.tool_columns}

    spec = resolve_key_spec(tool.tool_name)
    keyed: list[tuple[str, dict[str, str]]] = []
    if spec.kind == "synthetic_index":
        keyed = [(str(i + 1), r) for i, r in enumerate(rows)]
    else:
        keyed = [(format_composite_key(r, spec), r) for r in rows]

    seen_keys: set[str] = set()
    ordered_keyed: list[tuple[str, dict[str, str]]] = []
    for key, row in keyed:
        if key in seen_keys:
            continue
        seen_keys.add(key)
        ordered_keyed.append((key, row))

    if len(ordered_keyed) == 1:
        # Multiple raw rows, but all sharing the same key -- e.g. GeneProp
        # emits one row per matched PATHWAY STEP, so a single pathway with
        # several required steps produces several raw rows that all
        # collapse to one distinct key here. Render bare, same as the
        # len(rows) == 1 case above -- this is genuinely a single value,
        # not a multi-hit needing "key: val" wrapping. Without this check,
        # only key_cols would render bare (via the no-separator join
        # below) while every other column stayed wrapped as "key: val"
        # even though there was only one real entry.
        return {col: ordered_keyed[0][1].get(col, "") for col in tool.tool_columns}

    key_cols = set(spec.columns) if spec.kind in ("column",) else set()
    bare_extra_cols = BARE_VALUE_COLUMNS.get(tool.tool_name, set())

    merged: dict[str, str] = {}
    for col in tool.tool_columns:
        if col in key_cols:
            merged[col] = ";".join(key for key, _ in ordered_keyed if key)
        elif col in bare_extra_cols:
            # Bare deduplicated by the COLUMN's own values, not the primary key.
            seen_vals: set[str] = set()
            ordered_vals: list[str] = []
            for _, row in ordered_keyed:
                v = row.get(col, "")
                if v and v not in seen_vals:
                    seen_vals.add(v)
                    ordered_vals.append(v)
            merged[col] = ";".join(ordered_vals)
        elif tool.tool_name == "geneprop" and col == "GENEPROP_description":
            # Special render: "GenProp1235: Adenine and adenosine salvage
            # III(PARTIAL)" -- combines description + status into one
            # value per distinct GENEPROP_id (a generic per-column
            # id:value pairing would keep these separate).
            parts = []
            for key, row in ordered_keyed:
                desc = row.get(col, "")
                status = row.get("GENEPROP_status", "")
                parts.append(f"{key}: {desc}({status})" if key else f"{desc}({status})")
            merged[col] = "; ".join(parts)
        else:
            parts = []
            for key, row in ordered_keyed:
                val = row.get(col, "")
                parts.append(f"{key}: {val}" if key else val)
            merged[col] = "; ".join(parts)
    return merged


SHARED_INTERPRO_COLUMNS: tuple[str, ...] = (
    "INTERPRO_id", "INTERPRO_description", "INTERPRO_go_terms", "INTERPRO_pathways",
)


def aggregate_shared_interpro_columns(interpro_tables: list[ToolTable], fid: str) -> dict[str, str]:
    """INTERPRO_id/_description/_go_terms/_pathways are the cross-database
    InterPro entry reference -- the *same* literal column name appears in
    every interpro_<db> table (e.g. interpro_pfam_results.tsv and
    interpro_cdd_results.tsv both have their own "INTERPRO_id" column).
    Looping tool-by-tool and overwriting row[col] each time would silently
    keep only whichever interpro_<db> table happened to be processed last
    for this gene, discarding every other database's InterPro entries --
    a gene matching several distinct InterPro IDs across cdd/pfam/
    superfamily/etc. would otherwise lose all but the last one processed.
    This aggregates across ALL interpro_<db> tables' rows for this gene
    instead, deduped by InterPro ID, same "key: value" pairing convention
    as everywhere else once there's more than one distinct entry.
    """
    seen_ids: set[str] = set()
    ordered: list[dict[str, str]] = []
    for table in interpro_tables:
        for row in table.rows_by_feature.get(fid, []):
            ipr_id = row.get("INTERPRO_id", "")
            if not ipr_id or ipr_id in seen_ids:
                continue
            seen_ids.add(ipr_id)
            ordered.append(row)

    if not ordered:
        return {col: "" for col in SHARED_INTERPRO_COLUMNS}
    if len(ordered) == 1:
        return {col: ordered[0].get(col, "") for col in SHARED_INTERPRO_COLUMNS}

    merged: dict[str, str] = {}
    merged["INTERPRO_id"] = ";".join(r["INTERPRO_id"] for r in ordered)
    for col in ("INTERPRO_description", "INTERPRO_go_terms", "INTERPRO_pathways"):
        parts = [f"{r['INTERPRO_id']}: {r.get(col, '')}" for r in ordered]
        merged[col] = "; ".join(parts)
    return merged


def build_merged_rows(
    all_tables: list[ToolTable],
    organism_name_override: str,
    domain_override: str,
) -> tuple[list[dict], list[str]]:
    by_tool: dict[str, ToolTable] = {t.tool_name: t for t in all_tables}
    rast = by_tool.get("rasttk")

    all_fids: set[str] = set()
    for t in all_tables:
        all_fids.update(t.rows_by_feature)

    all_columns_order: list[str] = []
    seen_cols: set[str] = set()

    def register(col: str) -> None:
        if col not in seen_cols:
            seen_cols.add(col)
            all_columns_order.append(col)

    # gram_stain deliberately excluded -- ENVELOPE_envelope_type already
    # carries this (diderm-gram-negative-like / monoderm-gram-positive-
    # like / archaea), a separate gram_stain field would be redundant.
    for col in ("organism_name", "domain", "feature_id"):
        register(col)
    for col in RASTTK_IDENTITY_COLUMNS:
        register(col)
    for t in all_tables:
        for col in t.tool_columns:
            register(col)
    for col in ("ENVELOPE_envelope_type", "ENVELOPE_inference_basis", "ENVELOPE_evidence_json"):
        register(col)

    merged_rows: list[dict] = []
    for fid in sorted(all_fids):
        row: dict[str, str] = {
            "feature_id": fid,
            "organism_name": organism_name_override,
            "domain": domain_override,
        }
        for col in RASTTK_IDENTITY_COLUMNS:
            row[col] = ""
        if rast:
            rast_merged = merge_rows_for_gene(rast, fid)
            for col in rast.tool_columns:
                row[col] = rast_merged.get(col, "")
            # gene_id/gene_start/gene_end/na_length/aa_length/na_seq/aa_seq
            # are GLOBAL_COLUMNS -- deliberately excluded from tool_columns
            # (so they're not double-prefixed or multi-hit-wrapped), which
            # means rast_merged (built from tool_columns only) never has
            # them either. Pull them from rasttk's raw first row instead.
            rast_rows = rast.rows_by_feature.get(fid, [])
            if rast_rows:
                r0 = rast_rows[0]
                for col in RASTTK_IDENTITY_COLUMNS:
                    if r0.get(col):
                        row[col] = r0[col]
            if organism_name_override:
                row["organism_name"] = organism_name_override
            if domain_override:
                row["domain"] = domain_override

        for t in all_tables:
            if t.tool_name == "rasttk":
                continue
            tmerged = merge_rows_for_gene(t, fid)
            for col in t.tool_columns:
                row[col] = tmerged.get(col, "")

        interpro_tables = [t for t in all_tables if t.tool_name.startswith("interpro_")]
        if interpro_tables:
            row.update(aggregate_shared_interpro_columns(interpro_tables, fid))

        for col in ("ENVELOPE_envelope_type", "ENVELOPE_inference_basis", "ENVELOPE_evidence_json"):
            row.setdefault(col, "")
        for source_tool in ("deepsig", "psortb", "signalp4"):
            table = by_tool.get(source_tool)
            if not table:
                continue
            rows = table.rows_by_feature.get(fid, [])
            if not rows:
                continue
            r0 = rows[0]
            if r0.get("ENVELOPE_envelope_type"):
                for col in ("ENVELOPE_envelope_type", "ENVELOPE_inference_basis", "ENVELOPE_evidence_json"):
                    row[col] = r0.get(col, "")
                break

        merged_rows.append(row)

    return merged_rows, all_columns_order


# ─── Output ───────────────────────────────────────────────────────────────────

def write_merged_table(rows: list[dict], columns: list[str], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_manifest(all_tables: list[ToolTable], manifest_path: Path) -> None:
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    with open(manifest_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["tool_name", "source_file_path", "n_feature_ids", "n_tool_columns"])
        for t in all_tables:
            writer.writerow([t.tool_name, str(t.source_path), len(t.rows_by_feature), len(t.tool_columns)])


# ─── CLI ──────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input-root", required=True,
                   help="Genome output root, e.g. output/<genome>.")
    p.add_argument("--organism-name", default="")
    p.add_argument("--domain", default="", choices=["", "Bacteria", "Archaea", "Unknown"])
    p.add_argument("--exclude-tool", action="append", default=[])
    p.add_argument("--output", required=True, help="Path for the merged TSV.")
    p.add_argument("--manifest", default=None)
    return p.parse_args()


def main() -> None:
    args = parse_args()

    output_root = Path(args.input_root)
    if not output_root.is_dir():
        print(f"[merge-all-columns] ERROR: input root not a directory: {output_root}", file=sys.stderr)
        raise SystemExit(1)

    tool_table_pairs = discover_tool_tables(output_root, set(args.exclude_tool))
    if not tool_table_pairs:
        print("[merge-all-columns] ERROR: no per-tool tables found", file=sys.stderr)
        raise SystemExit(1)

    print(f"[merge-all-columns] {len(tool_table_pairs)} tool tables:")
    all_tables: list[ToolTable] = []
    for tool_name, src in tool_table_pairs:
        table = load_tool_table(tool_name, src)
        n_before = len(table.rows_by_feature)
        table = apply_tool_row_filter(table)
        all_tables.append(table)
        suffix = f" (filtered {n_before} -> {len(table.rows_by_feature)} feature_ids)" if n_before != len(table.rows_by_feature) else ""
        print(f"    {tool_name:14s} ← {src}  ({len(table.rows_by_feature)} feature_ids, {len(table.tool_columns)} cols){suffix}")

    merged_rows, columns = build_merged_rows(
        all_tables,
        organism_name_override=args.organism_name,
        domain_override=args.domain,
    )

    output_path = Path(args.output)
    write_merged_table(merged_rows, columns, output_path)
    print(f"[merge-all-columns] Wrote {len(merged_rows)} rows × {len(columns)} cols → {output_path}")

    if args.manifest:
        write_manifest(all_tables, Path(args.manifest))
        print(f"[merge-all-columns] Wrote manifest → {args.manifest}")

    print("[merge-all-columns] Done.")


if __name__ == "__main__":
    main()
