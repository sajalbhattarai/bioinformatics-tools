#!/usr/bin/env python3
"""update-fingerprint-database.py — margie_sb, incremental cross-genome
fingerprint-database update.

Takes ONE genome's labeled-genes-fingerprint-hash-label.tsv (phase12,
"pattern hash: <hash> || label: <canonical_label>" per gene) and merges it
into the SHARED, persistent fingerprint-database.tsv at
/depot/lindems/data/margie/fingerprint-database/ -- atomically, under an
exclusive file lock, so every genome's own phase12 run can update the same
shared file safely even if several genomes finish around the same time.

Mirrors build-here/.../phase10-fingerprinting/fingerprint/scripts/
process_fingerprint.py's own update_pangenome()/_do_update() pattern
exactly (fcntl.LOCK_EX on a sibling .lock file, write to .tmp, then
os.replace -- never partially-written, never two writers racing) -- that
script updates per-ORGANISM (one hash per organism, the set of its unique
canonical_labels); this one updates per-GENE (one hash per gene, the much
richer 15-tool evidence pattern from add-gene-fingerprint.py), but the
concurrency-safety shape is identical.

WHY incremental, not a full rebuild every time: re-reading every genome
processed so far just to add one more would get slower and slower as the
collection grows past dozens, then eventually hundreds of genomes. This
script only ever reads the existing pool once, merges in the one new
genome's contribution, and rewrites -- O(pool size + new genome's gene
count), not O(every genome ever processed).

Schema (fingerprint-database.tsv, one row per (pattern_id, label) pair,
same shape as the existing pangenome_fingerprints.tsv):
  pattern_id | fingerprint_hash | fingerprint_label
  | fingerprint_frequency | fingerprint_label_frequency | organisms

  fingerprint_frequency       : # genes (across every genome processed so
                                 far) that share this exact pattern
  fingerprint_label_frequency : # genes (across every pattern that decides
                                 on it) carrying this label
  organisms                   : pipe-separated organism names with at
                                 least one gene matching this pattern
"""
from __future__ import annotations

import argparse
import csv
import fcntl
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

csv.field_size_limit(10_000_000)

_HASH_LABEL_RE = re.compile(r"^pattern hash: (?P<hash>\S+) \|\| label: (?P<label>.*)$")


def parse_hash_label(fingerprint_value: str) -> tuple[str, str] | None:
    m = _HASH_LABEL_RE.match(fingerprint_value)
    if not m:
        return None
    return m.group("hash"), m.group("label")


def read_new_genome(hash_label_tsv: Path, organism: str) -> dict[str, int]:
    """Returns {hash: gene_count} for this one genome, plus stamps each
    hash's label as a side channel the caller merges separately (a single
    hash's label is assumed identical across every gene that has it,
    enforced by assign-canonical-label.py being a pure function of the
    same id/description fields the hash is computed from)."""
    counts: dict[str, int] = {}
    labels: dict[str, str] = {}
    n_unparsed = 0
    with open(hash_label_tsv, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            parsed = parse_hash_label(row.get("fingerprint", ""))
            if parsed is None:
                n_unparsed += 1
                continue
            h, label = parsed
            counts[h] = counts.get(h, 0) + 1
            labels[h] = label
    if n_unparsed:
        print(f"[update-fingerprint-database] WARNING: {n_unparsed} rows in {hash_label_tsv} "
              f"did not match the expected fingerprint format, skipped", file=sys.stderr)
    return counts, labels


def update_database(counts: dict[str, int], labels: dict[str, str], organism: str,
                     db_tsv: Path) -> None:
    db_tsv.parent.mkdir(parents=True, exist_ok=True)
    lock_path = db_tsv.with_suffix(".lock")
    with open(lock_path, "w") as lock_fh:
        fcntl.flock(lock_fh, fcntl.LOCK_EX)
        try:
            _do_update(counts, labels, organism, db_tsv)
        finally:
            fcntl.flock(lock_fh, fcntl.LOCK_UN)


def _do_update(new_counts: dict[str, int], new_labels: dict[str, str], organism: str,
                db_tsv: Path) -> None:
    """Inner update logic -- must be called under the exclusive lock above.

    patterns dict structure:
      hash -> {"pattern_id": str, "label": str, "gene_count": int, "organisms": set[str]}
    """
    patterns: dict[str, dict] = {}
    if db_tsv.exists() and db_tsv.stat().st_size > 0:
        with db_tsv.open(newline="", encoding="utf-8") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                h = row["fingerprint_hash"]
                patterns[h] = {
                    "pattern_id": row["pattern_id"],
                    "label": row["fingerprint_label"],
                    "gene_count": int(row["fingerprint_frequency"]),
                    "organisms": set(o for o in row.get("organisms", "").split("|") if o),
                }

    for h, count in new_counts.items():
        if h in patterns:
            patterns[h]["gene_count"] += count
            patterns[h]["organisms"].add(organism)
            # Defensive: identical hash should always mean identical label.
            if patterns[h]["label"] != new_labels[h]:
                print(f"[update-fingerprint-database] WARNING: hash {h} previously mapped to "
                      f"label '{patterns[h]['label']}', now seeing '{new_labels[h]}' from {organism} "
                      f"-- keeping the original, this shouldn't happen if assign-canonical-label.py "
                      f"is deterministic", file=sys.stderr)
        else:
            patterns[h] = {
                "pattern_id": f"pattern_{len(patterns) + 1}",
                "label": new_labels[h],
                "gene_count": count,
                "organisms": {organism},
            }

    label_freq: dict[str, int] = {}
    for pdata in patterns.values():
        label_freq[pdata["label"]] = label_freq.get(pdata["label"], 0) + pdata["gene_count"]

    tmp_path = db_tsv.with_suffix(".tmp")
    with tmp_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["pattern_id", "fingerprint_hash", "fingerprint_label",
                         "fingerprint_frequency", "fingerprint_label_frequency", "organisms"])
        for h, pdata in sorted(patterns.items(), key=lambda item: -item[1]["gene_count"]):
            writer.writerow([
                pdata["pattern_id"], h, pdata["label"], pdata["gene_count"],
                label_freq[pdata["label"]], "|".join(sorted(pdata["organisms"])),
            ])
    tmp_path.replace(db_tsv)

    meta_path = db_tsv.parent / "fingerprint-database-metadata.json"
    meta = {
        "total_organisms": len({o for p in patterns.values() for o in p["organisms"]}),
        "total_genes": sum(p["gene_count"] for p in patterns.values()),
        "total_patterns": len(patterns),
        "total_unique_labels": len(label_freq),
        "last_updated": datetime.now(timezone.utc).isoformat(),
        "last_updated_by": organism,
    }
    meta_path.write_text(json.dumps(meta, indent=2), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--hash-label-input", required=True,
                        help="this genome's labeled-genes-fingerprint-hash-label.tsv")
    parser.add_argument("--organism", required=True, help="organism identifier to stamp in the pool")
    parser.add_argument("--fingerprint-database", required=True,
                        help="shared fingerprint-database.tsv path (created if absent, updated in place)")
    args = parser.parse_args()

    hash_label_path = Path(args.hash_label_input)
    if not hash_label_path.is_file():
        print(f"[update-fingerprint-database] ERROR: input not found: {hash_label_path}", file=sys.stderr)
        raise SystemExit(1)

    db_path = Path(args.fingerprint_database)
    counts, labels = read_new_genome(hash_label_path, args.organism)
    print(f"[update-fingerprint-database] {args.organism}: {sum(counts.values())} genes, "
          f"{len(counts)} distinct patterns")
    print(f"[update-fingerprint-database] Updating: {db_path}")
    update_database(counts, labels, args.organism, db_path)
    print(f"[update-fingerprint-database] Done.")


if __name__ == "__main__":
    main()
