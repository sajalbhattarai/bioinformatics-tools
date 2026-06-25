#!/usr/bin/env python3
"""score-c3-pathway-coherence.py — margie_sb phase11 (scoring), metric
C3: GeneProp pathway coherence with operon-mates.

Reads consolidated-only-geneprop-pathway-status.tsv (phase9,
extract-geneprop-pathway-status.py's output) and
labeled-genes-operon-info.tsv (phase10, add-operon-info.py's output) --
BOTH READ-ONLY -- groups genes by operon_id, and checks whether each
gene shares a GeneProp pathway with at least one of its operon-mates.

Real corroborating evidence, not just a structural fact: a spurious
operon prediction is less likely to show functional coherence among its
members by chance, so finding a shared pathway is independent support
that BOTH the operon call and the underlying gene-function calls are
real. SEED-subsystem agreement (deferred -- no such column exists
anywhere in this pipeline currently) would be a second, complementary
half of this same idea, added later if that data source materialises.

Three real outcomes, not just match/no-match:
  0.75 -- this gene and >=1 operon-mate share a GeneProp pathway ID
          (checked, found, the corroborating case). Best-match
          aggregation, not unanimous -- only ONE corroborating
          operon-mate is needed; requiring every member of a large
          operon to share one pathway would almost never trigger even
          for genuine, well-annotated operons.
  0.0   -- this gene and/or its operon-mates DO carry GeneProp data, but
          none of it overlaps (checked, found NO overlap -- a real, if
          heuristic, negative signal).
  0.5   -- neutral: this gene is a singleton/non-coding (no operon-mates
          to check against at all), OR neither this gene nor any
          operon-mate carries any GeneProp data (nothing available to
          check) -- absence of a checkable comparison isn't evidence
          against the label, same "no_evidence is neutral, not a
          penalty" principle used in C4.

Output (labeled-genes-c3-pathway-coherence.tsv): identity columns,
c3_score, c3_shared_geneprop_ids (which pathway(s) corroborated, if any),
c3_formula (literal text explaining the outcome).
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

csv.field_size_limit(10_000_000)

_NEUTRAL = 0.5

_IDENTITY_COLUMNS = ["feature_id", "organism_name", "canonical_label", "label_source", "label_source_id"]


def parse_pathway_ids(geneprop_pathway_status: str) -> set[str]:
    """'GenProp0263:PARTIAL;GenProp1162:PARTIAL' -> {GenProp0263, GenProp1162}."""
    if not geneprop_pathway_status:
        return set()
    ids = set()
    for pair in geneprop_pathway_status.split(";"):
        if ":" in pair:
            ids.add(pair.split(":", 1)[0])
    return ids


def compute_c3(
    feature_id: str, operon_id: str, geneprop_by_gene: dict[str, set[str]],
    operon_members: dict[str, list[str]],
) -> tuple[float, str, str]:
    """Returns (score, shared_geneprop_ids, formula_text)."""
    if operon_id in ("NOT_IN_AN_OPERON", "NOT_APPLICABLE_NON_CODING"):
        return _NEUTRAL, "", f"no operon-mates to compare ({operon_id}) -> c3_score={_NEUTRAL:.4f} (neutral)"

    this_geneprop = geneprop_by_gene.get(feature_id, set())
    mates = [m for m in operon_members.get(operon_id, []) if m != feature_id]
    mate_geneprops = {m: geneprop_by_gene.get(m, set()) for m in mates}

    shared: set[str] = set()
    sharing_mates: list[str] = []
    for mate, mg in mate_geneprops.items():
        overlap = this_geneprop & mg
        if overlap:
            shared |= overlap
            sharing_mates.append(mate)

    if shared:
        formula = (f"shared GeneProp pathway(s) {sorted(shared)} with operon-mate(s) "
                  f"{sorted(sharing_mates)} -> c3_score=0.7500")
        return 0.75, ";".join(sorted(shared)), formula

    has_any_data = bool(this_geneprop) or any(mate_geneprops.values())
    if has_any_data:
        formula = (f"checked {len(mates)} operon-mate(s), no shared GeneProp pathway found "
                  f"-> c3_score=0.0000")
        return 0.0, "", formula

    formula = f"no GeneProp data on this gene or its {len(mates)} operon-mate(s) -> c3_score={_NEUTRAL:.4f} (neutral)"
    return _NEUTRAL, "", formula


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--geneprop-input", required=True,
                        help="extract-geneprop-pathway-status.py's output TSV "
                             "(consolidated-only-geneprop-pathway-status.tsv)")
    parser.add_argument("--operon-input", required=True, help="labeled-genes-operon-info.tsv")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    geneprop_path = Path(args.geneprop_input)
    operon_path = Path(args.operon_input)
    if not geneprop_path.is_file():
        print(f"[score-c3-pathway-coherence] ERROR: input not found: {geneprop_path}", file=sys.stderr)
        raise SystemExit(1)
    if not operon_path.is_file():
        print(f"[score-c3-pathway-coherence] ERROR: input not found: {operon_path}", file=sys.stderr)
        raise SystemExit(1)

    geneprop_by_gene: dict[str, set[str]] = {}
    with open(geneprop_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            fid = row.get("feature_id", "")
            if fid:
                geneprop_by_gene[fid] = parse_pathway_ids(row.get("geneprop_pathway_status", ""))

    operon_id_by_gene: dict[str, str] = {}
    operon_members: dict[str, list[str]] = {}
    with open(operon_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            fid = row.get("feature_id", "")
            if not fid:
                continue
            operon_id = row.get("operon_id", "")
            operon_id_by_gene[fid] = operon_id
            if operon_id.startswith("operon_"):
                operon_members.setdefault(operon_id, []).append(fid)

    out_columns = _IDENTITY_COLUMNS + ["c3_score", "c3_shared_geneprop_ids", "c3_formula"]

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    score_counts: dict[str, int] = {}
    n = 0
    with open(operon_path, newline="") as fh, open(output_path, "w", newline="") as out_fh:
        reader = csv.DictReader(fh, delimiter="\t")
        writer = csv.DictWriter(out_fh, fieldnames=out_columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in reader:
            fid = row.get("feature_id", "")
            operon_id = operon_id_by_gene.get(fid, "NOT_APPLICABLE_NON_CODING")
            c3, shared_ids, formula = compute_c3(fid, operon_id, geneprop_by_gene, operon_members)

            out_row = {col: row.get(col, "") for col in _IDENTITY_COLUMNS}
            out_row["c3_score"] = f"{c3:.4f}"
            out_row["c3_shared_geneprop_ids"] = shared_ids
            out_row["c3_formula"] = formula
            writer.writerow(out_row)
            score_counts[f"{c3:.4f}"] = score_counts.get(f"{c3:.4f}", 0) + 1
            n += 1

    print(f"[score-c3-pathway-coherence] Wrote {n} genes → {output_path}")
    for score, count in sorted(score_counts.items(), key=lambda kv: -kv[1]):
        print(f"    c3={score}  {count:6d} ({100.0 * count / n:.1f}%)")


if __name__ == "__main__":
    main()
