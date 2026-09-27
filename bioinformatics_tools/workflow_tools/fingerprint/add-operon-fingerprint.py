#!/usr/bin/env python3
"""add-operon-fingerprint.py — fingerprint stage, per-operon fingerprints.

Combines the gene-level hashes and labels of each operon's members into four
SHA-256 hashes: by evidence (gene pattern hashes) and by label, each in
genomic order and sorted. Evidence hashes are stricter; label hashes match
across species. Writes one row per gene, repeating its operon's fingerprint;
genes outside an operon get none.
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
    """Returns the first 16 hex characters of the string's SHA-256 (hashlib)."""
    return hashlib.sha256(s.encode("utf-8")).hexdigest()[:16]


def parse_hash_label(fingerprint_value: str) -> tuple[str, str] | None:
    """Extracts (hash, label) from "pattern hash: ... || label: ..." by regex, or None."""
    m = _HASH_LABEL_RE.match(fingerprint_value)
    if not m:
        return None
    return m.group("hash"), m.group("label")


def main() -> None:
    """Groups genes by operon with csv, hashes each operon's members and writes the per-gene table."""
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

    # operon_id -> [(position, feature_id)]
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

    # operon_id -> fingerprint string, "" when a member has no gene fingerprint
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
