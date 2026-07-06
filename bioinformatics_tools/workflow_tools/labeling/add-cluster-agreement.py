#!/usr/bin/env python3
"""add-cluster-agreement.py — margie_sb phase10 (labeling): tool-cluster
agreement.

Reads labeled-genes.tsv (assign-canonical-label.py's output, READ-ONLY)
and consolidated-merged-all-columns.tsv (for the underlying per-tool
accession columns), joins on feature_id, and writes its own standalone
table -- identity columns plus two distinct kinds of signal.

PURPOSE 1 -- collapsed cluster value, to fix scoring/score-c1-tool-
coverage.py's double-counting: PGAP, TIGRFAM, and NCBIFAM (via InterPro)
are three separate decision-tool slots in the trust hierarchy and in
C1's denominator, but they are not three independent sources -- PGAP's
own HMM library is built directly from TIGRFAM/NCBIfam's curated
models, and NCBIfam absorbed TIGRFAM outright. Confirmed empirically:
on a real genome, every gene where >=2 of these three had a hit, 100%
shared the same accession (modulo PGAP's own trailing ".N" version
suffix). tigrfam_cluster_* below collapses these three into one value/
status, so C1 can count this as ONE slot instead of three.

PURPOSE 2 -- confirmatory cross-reference, NOT consumed by C1 at all:
EGGNOG's own decision-tool vote (its best-hit description) is genuinely
independent -- a real DIAMOND/HMM search against eggNOG's own database.
But EGGNOG also reports a COG category (EGGNOG_COG_category, sometimes
a specific COG number embedded directly in EGGNOG_description) and a
KEGG KO number (EGGNOG_KEGG_ko) -- both LOOKED UP from its own
orthologous group's precomputed annotation table, not independently
searched. Comparing these against the standalone COG/KEGG tools' own
direct hits is real corroborating-or-contradicting evidence (same role
as the existing MEROPS/TCDB/DBCAN confirmatory check in assign-
canonical-label.py), but it must never be folded into C1's tool-
coverage count, since EGGNOG's own slot there is already counted via
its independent description.
Confirmed empirically these genuinely disagree often enough to matter:
COG vs EGGNOG's embedded COG# agree ~86% of the time (family-ID level,
when EGGNOG names one); KEGG vs EGGNOG_KEGG_ko agree ~80% of the time.
Neither is close to the TIGRFAM cluster's ~100%, since EGGNOG's value
is inherited at the whole-ortholog-group level, not the same per-gene
hit.

cdd_cog_overlap: checked directly against the same accession-prefix
question -- on every install checked so far, this InterPro build only
ever surfaces "cd#####"-prefixed CDD accessions, never the "COG####"-
formatted entries that also live inside NCBI's broader CDD/cddid.tbl
distribution. So this currently always reports False, but the check
itself stays general (a literal COG-prefixed token in
INTERPRO_CDD_id) in case a different InterPro/database version ever
surfaces one.
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

csv.field_size_limit(10_000_000)

_IDENTITY_COLUMNS = [
    "feature_id",
    "organism_name",
    "best_consensus_product_descriptor",
    "product_descriptor_source",
    "product_descriptor_source_id",
]

COG_NUMBER_RE = re.compile(r'\bCOG\d{4}\b')

TIGRFAM_CLUSTER_PREFERENCE = ["TIGRFAM", "PGAP", "NCBIFAM"]
TIGRFAM_CLUSTER_COLUMNS = {
    "TIGRFAM": "TIGRFAM_id",
    "PGAP": "PGAP_id",
    "NCBIFAM": "INTERPRO_NCBIFAM_id",
}


def normalize_accessions(raw: str) -> set[str]:
    """';'-joined accessions -> a set, stripping any trailing '.N' version
    suffix (PGAP_id carries one, e.g. 'TIGR02928.1'; the others don't)."""
    out: set[str] = set()
    for tok in raw.split(";"):
        tok = tok.strip()
        if not tok:
            continue
        head, _, tail = tok.rpartition(".")
        out.add(head if head and tail.isdigit() else tok)
    return out


def classify_agreement(evidence: dict[str, set[str]]) -> tuple[str, str, str]:
    """Generic per-tool value-set agreement classifier -- same shape as
    add-ec-consensus.py's classify_ec_agreement(), without that script's
    EC-specific wildcard-compatibility clustering (these accessions
    don't have a coarser/finer wildcard concept the way EC numbers do,
    so plain set intersection/union is the right level of complexity
    here). Returns (status, consensus_value, supporting_tools)."""
    if not evidence:
        return "no_evidence", "", ""
    tools = list(evidence.keys())
    if len(tools) == 1:
        only = tools[0]
        return "single_source", ";".join(sorted(evidence[only])), only
    sets = list(evidence.values())
    intersection = set.intersection(*sets)
    union = set.union(*sets)
    if intersection and intersection == union:
        return "full_consensus", ";".join(sorted(intersection)), ";".join(sorted(tools))
    if intersection:
        supporting = [t for t in tools if intersection <= evidence[t]]
        return "majority_consensus", ";".join(sorted(intersection)), ";".join(sorted(supporting))
    return "conflicting", "", ";".join(sorted(tools))


def compute_tigrfam_cluster(merged_row: dict[str, str]) -> tuple[str, str, str, str]:
    """Returns (cluster_value, cluster_source, cluster_status, formula)."""
    evidence: dict[str, set[str]] = {}
    for tool in TIGRFAM_CLUSTER_PREFERENCE:
        vals = normalize_accessions(merged_row.get(TIGRFAM_CLUSTER_COLUMNS[tool], ""))
        if vals:
            evidence[tool] = vals

    status, _, _ = classify_agreement(evidence)
    if status == "no_evidence":
        return "", "", status, "no TIGRFAM/PGAP/NCBIfam(InterPro) hit"

    source = next(t for t in TIGRFAM_CLUSTER_PREFERENCE if t in evidence)
    value = sorted(evidence[source])[0]
    parts = [f"{t}={'|'.join(sorted(v))}" for t, v in evidence.items()]
    formula = f"{'; '.join(parts)} -> {status} (representative: {source}={value})"
    return value, source, status, formula


def compute_cog_crossref(merged_row: dict[str, str]) -> tuple[str, str, str, str]:
    """Returns (crossref_value, crossref_level, status, formula).
    crossref_level is 'family_id' when EGGNOG_description names a
    specific COG number, else 'category_letter' when only the coarser
    EGGNOG_COG_category is available."""
    cog_id = merged_row.get("COG_id", "")
    cog_letter = merged_row.get("COG_func_letter", "")
    eggnog_desc = merged_row.get("EGGNOG_description", "")
    eggnog_cogs = set(COG_NUMBER_RE.findall(eggnog_desc))
    eggnog_letter = merged_row.get("EGGNOG_COG_category", "")

    if eggnog_cogs:
        evidence: dict[str, set[str]] = {"EGGNOG": eggnog_cogs}
        if cog_id:
            evidence["COG"] = set(cog_id.split(";"))
        status, _, _ = classify_agreement(evidence)
        formula = (f"COG_id={cog_id or '(none)'}; EGGNOG_description cites {sorted(eggnog_cogs)} "
                  f"-> {status} (family-ID level)")
        return ";".join(sorted(eggnog_cogs)), "family_id", status, formula

    if eggnog_letter and eggnog_letter != "-":
        evidence = {"EGGNOG": set(eggnog_letter)}
        if cog_letter:
            evidence["COG"] = set(cog_letter)
        status, _, _ = classify_agreement(evidence)
        formula = (f"COG_func_letter={cog_letter or '(none)'}; EGGNOG_COG_category={eggnog_letter} "
                  f"-> {status} (category-letter level, no specific COG# named)")
        return eggnog_letter, "category_letter", status, formula

    return "", "", "no_evidence", "EGGNOG gave no COG information at all"


def compute_kegg_crossref(merged_row: dict[str, str]) -> tuple[str, str, str]:
    """Returns (crossref_value, status, formula)."""
    kegg_id = merged_row.get("KEGG_id", "")
    eggnog_ko_raw = merged_row.get("EGGNOG_KEGG_ko", "")
    eggnog_ko = {t.replace("ko:", "").strip() for t in eggnog_ko_raw.split(",")
                if t.strip() and t.strip() != "-"}

    evidence: dict[str, set[str]] = {}
    if kegg_id:
        evidence["KEGG"] = set(kegg_id.split(";"))
    if eggnog_ko:
        evidence["EGGNOG"] = eggnog_ko
    status, _, _ = classify_agreement(evidence)
    formula = f"KEGG_id={kegg_id or '(none)'}; EGGNOG_KEGG_ko={eggnog_ko_raw or '(none)'} -> {status}"
    return ";".join(sorted(eggnog_ko)), status, formula


def compute_cdd_cog_overlap(merged_row: dict[str, str]) -> str:
    """True only if INTERPRO_CDD_id itself carries a COG-formatted token --
    on every install checked so far it doesn't (this InterPro build only
    ever surfaces 'cd#####'-prefixed CDD accessions), so this currently
    always reports False, but the check stays general."""
    cdd_id = merged_row.get("INTERPRO_CDD_id", "")
    return "True" if any(tok.strip().startswith("COG") for tok in cdd_id.split(";")) else "False"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--labeled-input", required=True, help="labeled-genes.tsv")
    parser.add_argument("--merged-input", required=True,
                        help="consolidated-merged-all-columns.tsv (for the underlying per-tool ID columns)")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    labeled_path = Path(args.labeled_input)
    merged_path = Path(args.merged_input)
    if not labeled_path.is_file():
        print(f"[add-cluster-agreement] ERROR: input not found: {labeled_path}", file=sys.stderr)
        raise SystemExit(1)
    if not merged_path.is_file():
        print(f"[add-cluster-agreement] ERROR: input not found: {merged_path}", file=sys.stderr)
        raise SystemExit(1)

    merged_by_gene: dict[str, dict[str, str]] = {}
    with open(merged_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            fid = row.get("feature_id", "")
            if fid:
                merged_by_gene[fid] = row

    out_columns = _IDENTITY_COLUMNS + [
        "tigrfam_cluster_value", "tigrfam_cluster_source", "tigrfam_cluster_status", "tigrfam_cluster_formula",
        "cog_crossref_value", "cog_crossref_level", "cog_crossref_status", "cog_crossref_formula",
        "kegg_crossref_value", "kegg_crossref_status", "kegg_crossref_formula",
        "cdd_cog_overlap",
    ]

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    tigrfam_status_counts: dict[str, int] = {}
    cog_status_counts: dict[str, int] = {}
    kegg_status_counts: dict[str, int] = {}
    n = 0
    with open(labeled_path, newline="") as fh, open(output_path, "w", newline="") as out_fh:
        reader = csv.DictReader(fh, delimiter="\t")
        writer = csv.DictWriter(out_fh, fieldnames=out_columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in reader:
            fid = row.get("feature_id", "")
            merged_row = merged_by_gene.get(fid, {})

            tigrfam_value, tigrfam_source, tigrfam_status, tigrfam_formula = compute_tigrfam_cluster(merged_row)
            cog_value, cog_level, cog_status, cog_formula = compute_cog_crossref(merged_row)
            kegg_value, kegg_status, kegg_formula = compute_kegg_crossref(merged_row)
            cdd_overlap = compute_cdd_cog_overlap(merged_row)

            out_row = {col: row.get(col, "") for col in _IDENTITY_COLUMNS}
            out_row["tigrfam_cluster_value"] = tigrfam_value
            out_row["tigrfam_cluster_source"] = tigrfam_source
            out_row["tigrfam_cluster_status"] = tigrfam_status
            out_row["tigrfam_cluster_formula"] = tigrfam_formula
            out_row["cog_crossref_value"] = cog_value
            out_row["cog_crossref_level"] = cog_level
            out_row["cog_crossref_status"] = cog_status
            out_row["cog_crossref_formula"] = cog_formula
            out_row["kegg_crossref_value"] = kegg_value
            out_row["kegg_crossref_status"] = kegg_status
            out_row["kegg_crossref_formula"] = kegg_formula
            out_row["cdd_cog_overlap"] = cdd_overlap
            writer.writerow(out_row)

            tigrfam_status_counts[tigrfam_status] = tigrfam_status_counts.get(tigrfam_status, 0) + 1
            cog_status_counts[cog_status] = cog_status_counts.get(cog_status, 0) + 1
            kegg_status_counts[kegg_status] = kegg_status_counts.get(kegg_status, 0) + 1
            n += 1

    print(f"[add-cluster-agreement] Wrote {n} genes -> {output_path}")
    print("  tigrfam_cluster_status:")
    for status, count in sorted(tigrfam_status_counts.items(), key=lambda kv: -kv[1]):
        print(f"    {status:20s} {count:6d} ({100.0 * count / n:.1f}%)")
    print("  cog_crossref_status:")
    for status, count in sorted(cog_status_counts.items(), key=lambda kv: -kv[1]):
        print(f"    {status:20s} {count:6d} ({100.0 * count / n:.1f}%)")
    print("  kegg_crossref_status:")
    for status, count in sorted(kegg_status_counts.items(), key=lambda kv: -kv[1]):
        print(f"    {status:20s} {count:6d} ({100.0 * count / n:.1f}%)")


if __name__ == "__main__":
    main()
