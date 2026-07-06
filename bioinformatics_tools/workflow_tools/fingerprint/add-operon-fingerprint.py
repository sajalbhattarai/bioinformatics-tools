#!/usr/bin/env python3
"""add-operon-fingerprint.py — margie_sb phase12 (fingerprint), per-operon
fingerprint.

Reads labeled-genes-operon-info.tsv (phase10) and labeled-genes-fingerprint-
hash-label.tsv (phase12, this same genome's own gene-level fingerprints) --
both READ-ONLY -- and composes each operon's member genes' own pattern
hashes/labels into FOUR operon-level signals, kept separate because each
answers a different question:

  by EVIDENCE (gene_pattern_hash):
    ordered hash      hash of [gene1_hash, gene2_hash, ...] in the genes'
                       actual genomic order (operon_gene_position_in_operon).
                       "has this exact arrangement -- same genes, same
                       order, same underlying tool evidence -- been seen
                       before?" Strict: two homologous operons in different
                       species will usually fail this even when the
                       biology is identical, since the specific tool hit
                       behind each gene's label can differ slightly
                       species to species.
    composition hash  hash of the same set, sorted. Same evidence-level
                       strictness, order-independent.

  by LABEL (canonical_label):
    ordered hash      hash of [label1, label2, ...] in genomic order.
                       "has this same DECIDED FUNCTION sequence been seen
                       before?" -- the one that actually generalizes across
                       species, since it tolerates different tools/specific
                       hits landing on the same final label.
    composition hash  hash of the same labels, sorted.

Since gene-level hash -> label is deterministic in practice (assign-
canonical-label.py is a pure function of the same id/description fields
the hash is computed from), the evidence-based hashes are always at least
as strict as the label-based ones -- they can only narrow what counts as
a match, never widen it.

Singletons (operon_id == NOT_IN_AN_OPERON) get no operon fingerprint at
all -- there's nothing to compose for a lone gene, same "blank, not a
penalty" principle used for C2/C3 throughout phase11.

Output (labeled-genes-operon-fingerprint.tsv): one row PER GENE, not per
operon -- every member of the same operon carries the identical operon
fingerprint value, repeated. Matches how operon_id/operon_probability
already attach per-gene in labeled-genes-operon-info.tsv rather than as a
separate operon-indexed table, so a consumer can look up "this gene's
operon fingerprint" directly by feature_id with no second join.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import re
import sys
from pathlib import Path

csv.field_size_limit(10_000_000)

_HASH_LABEL_RE = re.compile(r"^pattern hash: (?P<hash>\S+) \|\| label: (?P<label>.*)$")
_NOT_IN_OPERON = "NOT_IN_AN_OPERON"


def _hash16(s: str) -> str:
    return hashlib.sha256(s.encode("utf-8")).hexdigest()[:16]


def parse_hash_label(fingerprint_value: str) -> tuple[str, str] | None:
    m = _HASH_LABEL_RE.match(fingerprint_value)
    if not m:
        return None
    return m.group("hash"), m.group("label")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--operon-input", required=True, help="labeled-genes-operon-info.tsv")
    parser.add_argument("--hash-label-input", required=True,
                        help="this genome's own labeled-genes-fingerprint-hash-label.tsv")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    operon_path = Path(args.operon_input)
    hash_label_path = Path(args.hash_label_input)
    for p in (operon_path, hash_label_path):
        if not p.is_file():
            print(f"[add-operon-fingerprint] ERROR: input not found: {p}", file=sys.stderr)
            raise SystemExit(1)

    gene_hash_label: dict[str, tuple[str, str]] = {}
    with open(hash_label_path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            parsed = parse_hash_label(row.get("fingerprint", ""))
            if parsed:
                gene_hash_label[row["feature_id"]] = parsed

    # operon_id -> list of (position, feature_id), built up while reading
    operon_members: dict[str, list[tuple[int, str]]] = {}
    operon_id_by_gene: dict[str, str] = {}
    rows: list[dict[str, str]] = []
    with open(operon_path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            rows.append(row)
            fid = row["feature_id"]
            oid = row.get("operon_id", _NOT_IN_OPERON)
            operon_id_by_gene[fid] = oid
            if oid != _NOT_IN_OPERON:
                try:
                    pos = int(row.get("operon_gene_position_in_operon", "0"))
                except ValueError:
                    pos = 0
                operon_members.setdefault(oid, []).append((pos, fid))

    # operon_id -> formatted fingerprint string (or "" if any member lacks a
    # gene-level fingerprint, which shouldn't happen but isn't load-bearing
    # enough to crash over)
    operon_fingerprint: dict[str, str] = {}
    for oid, members in operon_members.items():
        members.sort(key=lambda t: t[0])
        ordered_fids = [fid for _, fid in members]
        if not all(fid in gene_hash_label for fid in ordered_fids):
            operon_fingerprint[oid] = ""
            continue
        ordered_hashes = [gene_hash_label[fid][0] for fid in ordered_fids]
        ordered_labels = [gene_hash_label[fid][1] for fid in ordered_fids]

        evidence_ordered_hash = _hash16(" | ".join(ordered_hashes))
        evidence_composition_hash = _hash16(" | ".join(sorted(ordered_hashes)))
        label_ordered_hash = _hash16(" | ".join(ordered_labels))
        label_composition_hash = _hash16(" | ".join(sorted(ordered_labels)))
        members_in_order = " -> ".join(ordered_labels)
        gene_pattern_hashes = " | ".join(ordered_hashes)

        operon_fingerprint[oid] = (
            f"operon hash by evidence (ordered): {evidence_ordered_hash} || "
            f"operon hash by evidence (composition): {evidence_composition_hash} || "
            f"operon hash by label (ordered): {label_ordered_hash} || "
            f"operon hash by label (composition): {label_composition_hash} || "
            f"members (in order): {members_in_order} || "
            f"gene_pattern_hashes: {gene_pattern_hashes}"
        )

    out_columns = ["organism_name", "feature_id", "operon_id", "operon_fingerprint"]
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    n = 0
    n_operonic = 0
    with open(output_path, "w", newline="") as out_fh:
        writer = csv.DictWriter(out_fh, fieldnames=out_columns, delimiter="\t")
        writer.writeheader()
        for row in rows:
            fid = row["feature_id"]
            oid = operon_id_by_gene.get(fid, _NOT_IN_OPERON)
            writer.writerow({
                "organism_name": row.get("organism_name", ""),
                "feature_id": fid,
                "operon_id": oid,
                "operon_fingerprint": operon_fingerprint.get(oid, ""),
            })
            n += 1
            if oid != _NOT_IN_OPERON:
                n_operonic += 1

    print(f"[add-operon-fingerprint] Wrote {n} genes → {output_path}")
    print(f"    in a real operon: {n_operonic} ({100.0*n_operonic/n:.1f}%), "
          f"across {len(operon_members)} distinct operons")


if __name__ == "__main__":
    main()
