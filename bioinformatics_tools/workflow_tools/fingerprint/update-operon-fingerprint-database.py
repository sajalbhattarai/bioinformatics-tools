#!/usr/bin/env python3
"""update-operon-fingerprint-database.py — merges one genome's operon fingerprints into the shared pools.

Counts each distinct operon once and adds it to four pools, keyed by the
evidence/label ordered/composition hashes from add-operon-fingerprint.py,
using the same fcntl lock and .tmp-then-rename write as
update-fingerprint-database.py. Columns: pattern_id, fingerprint_hash,
members_in_order, fingerprint_frequency (operons with this hash),
fingerprint_label_frequency (operons with this member sequence), organisms.
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

_OPERON_FP_RE = re.compile(
    r"^operon hash by evidence \(ordered\): (?P<evidence_ordered>\S+) \|\| "
    r"operon hash by evidence \(composition\): (?P<evidence_composition>\S+) \|\| "
    r"operon hash by label \(ordered\): (?P<label_ordered>\S+) \|\| "
    r"operon hash by label \(composition\): (?P<label_composition>\S+) \|\| "
    r"members \(in order\): (?P<members>.*?) \|\| "
    r"gene_pattern_hashes: .*$"
)


def parse_operon_fingerprint(value: str) -> dict | None:
    """Parses an operon fingerprint string by regex into its four hashes and member list, or None."""
    m = _OPERON_FP_RE.match(value)
    if not m:
        return None
    return m.groupdict()


def read_new_genome(operon_fp_tsv: Path) -> dict[str, dict]:
    """Reads the per-gene operon fingerprint TSV with csv and returns one parsed fingerprint per operon_id."""
    operons: dict[str, dict] = {}
    n_unparsed = 0
    with open(operon_fp_tsv, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            oid = row.get("operon_id", "")
            if oid == "NOT_IN_AN_OPERON" or oid in operons:
                continue
            parsed = parse_operon_fingerprint(row.get("operon_fingerprint", ""))
            if parsed is None:
                if row.get("operon_fingerprint", ""):
                    n_unparsed += 1
                continue
            operons[oid] = parsed
    if n_unparsed:
        print(f"[update-operon-fingerprint-database] WARNING: {n_unparsed} non-empty operon_fingerprint "
              f"values did not match the expected format, skipped", file=sys.stderr)
    return operons


def update_pool(operons: dict[str, dict], organism: str, db_tsv: Path, key: str) -> None:
    """Updates one pool keyed on the `key` hash, holding an exclusive fcntl lock."""
    db_tsv.parent.mkdir(parents=True, exist_ok=True)
    lock_path = db_tsv.with_suffix(".lock")
    with open(lock_path, "w") as lock_fh:
        fcntl.flock(lock_fh, fcntl.LOCK_EX)
        try:
            _do_update(operons, organism, db_tsv, key)
        finally:
            fcntl.flock(lock_fh, fcntl.LOCK_UN)


def _do_update(new_operons: dict[str, dict], organism: str, db_tsv: Path, key: str) -> None:
    """Merges operon counts into a pool TSV, rewrites it via .tmp and rename, and writes metadata JSON."""
    patterns: dict[str, dict] = {}
    if db_tsv.exists() and db_tsv.stat().st_size > 0:
        with db_tsv.open(newline="", encoding="utf-8") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                h = row["fingerprint_hash"]
                patterns[h] = {
                    "pattern_id": row["pattern_id"],
                    "members_in_order": row["members_in_order"],
                    "operon_count": int(row["fingerprint_frequency"]),
                    "organisms": set(o for o in row.get("organisms", "").split("|") if o),
                }

    # Each distinct operon counts once, even when several in one genome share a hash.
    for oid, parsed in new_operons.items():
        h = parsed[key]
        members = parsed["members"]
        if h in patterns:
            patterns[h]["operon_count"] += 1
            patterns[h]["organisms"].add(organism)
        else:
            patterns[h] = {
                "pattern_id": f"pattern_{len(patterns) + 1}",
                "members_in_order": members,
                "operon_count": 1,
                "organisms": {organism},
            }

    label_freq: dict[str, int] = {}
    for pdata in patterns.values():
        label_freq[pdata["members_in_order"]] = label_freq.get(pdata["members_in_order"], 0) + pdata["operon_count"]

    tmp_path = db_tsv.with_suffix(".tmp")
    with tmp_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["pattern_id", "fingerprint_hash", "members_in_order",
                         "fingerprint_frequency", "fingerprint_label_frequency", "organisms"])
        for h, pdata in sorted(patterns.items(), key=lambda item: -item[1]["operon_count"]):
            writer.writerow([
                pdata["pattern_id"], h, pdata["members_in_order"], pdata["operon_count"],
                label_freq[pdata["members_in_order"]], "|".join(sorted(pdata["organisms"])),
            ])
    tmp_path.replace(db_tsv)

    meta_path = db_tsv.parent / f"operon-fingerprint-database-{key}-metadata.json"
    meta = {
        "total_organisms": len({o for p in patterns.values() for o in p["organisms"]}),
        "total_operons": sum(p["operon_count"] for p in patterns.values()),
        "total_patterns": len(patterns),
        "total_unique_member_sequences": len(label_freq),
        "last_updated": datetime.now(timezone.utc).isoformat(),
        "last_updated_by": organism,
    }
    meta_path.write_text(json.dumps(meta, indent=2), encoding="utf-8")


_POOLS = [
    ("evidence_ordered", "--evidence-ordered-database"),
    ("evidence_composition", "--evidence-composition-database"),
    ("label_ordered", "--label-ordered-database"),
    ("label_composition", "--label-composition-database"),
]


def main() -> None:
    """Reads the genome's operon fingerprints and updates all four pools."""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--operon-fingerprint-input", required=True,
                        help="this genome's labeled-genes-operon-fingerprint.tsv")
    parser.add_argument("--organism", required=True, help="organism identifier to stamp in the pools")
    parser.add_argument("--evidence-ordered-database", required=True)
    parser.add_argument("--evidence-composition-database", required=True)
    parser.add_argument("--label-ordered-database", required=True)
    parser.add_argument("--label-composition-database", required=True)
    args = parser.parse_args()

    operon_fp_path = Path(args.operon_fingerprint_input)
    if not operon_fp_path.is_file():
        print(f"[update-operon-fingerprint-database] ERROR: input not found: {operon_fp_path}", file=sys.stderr)
        raise SystemExit(1)

    operons = read_new_genome(operon_fp_path)
    print(f"[update-operon-fingerprint-database] {args.organism}: {len(operons)} distinct operons")

    for key, flag in _POOLS:
        db_path = getattr(args, flag.lstrip("-").replace("-", "_"))
        update_pool(operons, args.organism, Path(db_path), key)
        print(f"[update-operon-fingerprint-database] Updated {key} pool: {db_path}")
    print(f"[update-operon-fingerprint-database] Done.")


if __name__ == "__main__":
    main()
