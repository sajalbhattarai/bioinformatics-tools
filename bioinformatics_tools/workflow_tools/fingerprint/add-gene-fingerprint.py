#!/usr/bin/env python3
"""add-gene-fingerprint.py — margie_sb phase12 (fingerprint), per-gene fingerprint.

Runs after scoring (phase11), not after labeling -- the full-with-scores
output needs labeled-genes-confidence-final.tsv, so fingerprinting can't
start until scoring has finished.

Reads labeled-genes.tsv (phase10) and labeled-genes-confidence-final.tsv
(phase11/scoring) -- both READ-ONLY -- and writes FIVE slimmed derived
views, one per useful combination of the three core pieces of information
a gene's fingerprint carries: its raw evidence pattern, that pattern's
hash, and the label decided from it. None of the four carry the ~85 wide
per-tool evidence columns labeled-genes.tsv has -- those already got
distilled into the fingerprint values themselves.

THE THREE CORE PIECES:
  hash         SHA-256 of the raw fingerprint values, truncated to 16 hex
               chars (same truncation the existing per-organism
               fingerprint.sif container already uses).
  label        canonical_label, exactly as labeling decided it.
  fingerprint  the raw values themselves: a FIXED-POSITION, pipe-joined
               list -- RAST_description first, then {tool}_all_ids,
               {tool}_all_descriptions for each decision tool in a fixed
               order (PGAP, TIGRFAM, HAMAP, NCBIFAM, PIRSF, UNIPROT, PFAM,
               CDD, KEGG, EGGNOG, COG, MEROPS, TCDB, DBCAN). Slots stay
               empty rather than being skipped when a tool has no hit, so
               position is comparable across genes regardless of which
               tools fired -- two genes are only an exact fingerprint
               match if every tool that fired (and every tool that didn't)
               lines up.

FOUR OUTPUT FILES, each combining two-or-three of those pieces (every
combination that includes at least two -- a bare hash, bare label, or bare
fingerprint alone isn't useful on its own):

  labeled-genes-fingerprint-hash-pattern.tsv
      "pattern hash: <hash> || fingerprint: <values>"
      Clusters genes by identical raw evidence regardless of what label
      won -- doesn't care what got decided, only what every tool actually
      found.

  labeled-genes-fingerprint-hash-label.tsv
      "pattern hash: <hash> || label: <canonical_label>"
      Compact hash-to-label lookup with no raw values repeated -- the
      shape a future cross-genome fingerprint-database dedup table wants.

  labeled-genes-fingerprint-label-pattern.tsv
      "label: <canonical_label> || fingerprint: <values>"
      Human-readable audit view -- label and literal evidence side by
      side, no hash to look up separately.

  labeled-genes-fingerprint-full.tsv
      "pattern hash: <hash> || label: <canonical_label> || fingerprint: <values>"
      The complete, self-contained record.

A FIFTH FILE adds the C1-C4/confidence_score layer on top of the full
record:

  labeled-genes-fingerprint-full-with-scores.tsv
      "pattern hash: <hash> || label: <canonical_label> || scores: C1:..
      |C2:.. |C3:.. |C4:.. |CONFIDENCE:<score>:<tier> || fingerprint: <values>"
      Same as full, plus the gene's confidence-score breakdown -- needs
      labeled-genes-confidence-final.tsv (phase11/scoring), so non-coding
      features (no scoring row at all) get "scores: " left empty rather
      than fabricating a score that was never computed.

Non-coding features still get a fingerprint (RAST_description alone, if
present) in the other four files -- fingerprinting only depends on
labeled-genes.tsv, not on scoring, so there's no reason to blank those out
the way the scores layer is blanked for non-coding rows.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import sys
from pathlib import Path

csv.field_size_limit(10_000_000)

_KEPT_COLUMNS = [
    "organism_name", "feature_id", "domain", "gene_id", "gene_start", "gene_end",
    "RAST_feature_type", "RAST_strand",
]

# Fixed positional order -- RAST has no id, so it contributes one slot
# (its description); every other tool contributes two slots (id, then
# description), always present even when empty, so two genes' fingerprint
# values line up slot-for-slot regardless of which tools happened to fire.
# Same tool list/order assign-canonical-label.py's trust hierarchy walks.
_ALL_IDS_DESC_TOOLS = [
    "PGAP", "TIGRFAM", "HAMAP", "NCBIFAM", "PIRSF", "UNIPROT", "PFAM",
    "CDD", "KEGG", "EGGNOG", "COG", "MEROPS", "TCDB", "DBCAN",
]


def _hash16(s: str) -> str:
    return hashlib.sha256(s.encode("utf-8")).hexdigest()[:16]


def _all_slots_empty(fingerprint_values: str) -> bool:
    """True if every "{field}: {value}" slot's value half is empty --
    i.e. no tool fired at all, not even RAST's bare gene-calling
    description."""
    return all(slot.split(": ", 1)[-1] == "" for slot in fingerprint_values.split(" | "))


def build_fingerprint_values(row: dict[str, str]) -> str:
    """Each slot is "{field_name}: {value}", not a bare value -- field
    name always present even when the value is empty, so a slot is
    self-describing on its own, not just by position."""
    slots = [f"RAST_description: {row.get('RAST_description', '').strip()}"]
    for tool in _ALL_IDS_DESC_TOOLS:
        slots.append(f"{tool}_id: {row.get(f'{tool}_all_ids', '').strip()}")
        slots.append(f"{tool}_description: {row.get(f'{tool}_all_descriptions', '').strip()}")
    return " | ".join(slots)


def build_scores_token(score_row: dict[str, str] | None) -> str:
    if score_row is None:
        return ""
    return (
        f"C1:{score_row.get('c1_score', '')}"
        f"|C2:{score_row.get('c2_score_from_operon_probability', '')}"
        f"|C3:{score_row.get('c3_score', '')}"
        f"|C4:{score_row.get('c4_score', '')}"
        f"|CONFIDENCE:{score_row.get('confidence_score', '')}:{score_row.get('confidence_score_tier', '')}"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--labeled-input", required=True, help="labeled-genes.tsv")
    parser.add_argument("--confidence-final-input", required=True,
                        help="labeled-genes-confidence-final.tsv (phase11/scoring) -- "
                             "only consumed by --output-full-with-scores")
    parser.add_argument("--output-hash-pattern", required=True)
    parser.add_argument("--output-hash-label", required=True)
    parser.add_argument("--output-label-pattern", required=True)
    parser.add_argument("--output-full", required=True)
    parser.add_argument("--output-full-with-scores", required=True)
    args = parser.parse_args()

    labeled_path = Path(args.labeled_input)
    confidence_path = Path(args.confidence_final_input)
    for p in (labeled_path, confidence_path):
        if not p.is_file():
            print(f"[add-gene-fingerprint] ERROR: input not found: {p}", file=sys.stderr)
            raise SystemExit(1)

    score_by_fid: dict[str, dict[str, str]] = {}
    with open(confidence_path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            fid = row.get("feature_id", "")
            if fid:
                score_by_fid[fid] = row

    outputs = {
        "hash_pattern": Path(args.output_hash_pattern),
        "hash_label": Path(args.output_hash_label),
        "label_pattern": Path(args.output_label_pattern),
        "full": Path(args.output_full),
        "full_with_scores": Path(args.output_full_with_scores),
    }
    for p in outputs.values():
        p.parent.mkdir(parents=True, exist_ok=True)

    out_columns = _KEPT_COLUMNS + ["fingerprint"]

    n = 0
    no_hit_count = 0
    with open(labeled_path, newline="") as fh, \
         open(outputs["hash_pattern"], "w", newline="") as f_hp, \
         open(outputs["hash_label"], "w", newline="") as f_hl, \
         open(outputs["label_pattern"], "w", newline="") as f_lp, \
         open(outputs["full"], "w", newline="") as f_full, \
         open(outputs["full_with_scores"], "w", newline="") as f_fws:

        reader = csv.DictReader(fh, delimiter="\t")
        writers = {
            key: csv.DictWriter(f, fieldnames=out_columns, delimiter="\t", extrasaction="ignore")
            for key, f in (("hash_pattern", f_hp), ("hash_label", f_hl),
                           ("label_pattern", f_lp), ("full", f_full),
                           ("full_with_scores", f_fws))
        }
        for w in writers.values():
            w.writeheader()

        for row in reader:
            values = build_fingerprint_values(row)
            h = _hash16(values)
            label = row.get("best_consensus_product_descriptor", "")
            identity = {col: row.get(col, "") for col in _KEPT_COLUMNS}
            scores = build_scores_token(score_by_fid.get(row.get("feature_id", "")))

            writers["hash_pattern"].writerow({**identity, "fingerprint": f"pattern hash: {h} || fingerprint: {values}"})
            writers["hash_label"].writerow({**identity, "fingerprint": f"pattern hash: {h} || label: {label}"})
            writers["label_pattern"].writerow({**identity, "fingerprint": f"label: {label} || fingerprint: {values}"})
            writers["full"].writerow({**identity, "fingerprint": f"pattern hash: {h} || label: {label} || fingerprint: {values}"})
            full_with_scores = (f"pattern hash: {h} || label: {label} || scores: {scores} || fingerprint: {values}"
                                 if scores else f"pattern hash: {h} || label: {label} || fingerprint: {values}")
            writers["full_with_scores"].writerow({**identity, "fingerprint": full_with_scores})

            if _all_slots_empty(values):
                no_hit_count += 1
            n += 1

    print(f"[add-gene-fingerprint] Wrote {n} genes x 5 files:")
    for key, p in outputs.items():
        print(f"    {key}: {p}")
    print(f"    genes with zero tool hits (empty fingerprint values): {no_hit_count} ({100.0*no_hit_count/n:.1f}%)")


if __name__ == "__main__":
    main()
