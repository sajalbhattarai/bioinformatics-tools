#!/usr/bin/env python3
"""Update the depot-hosted OCC reference with one organism (idempotent).

This runs before per-organism C3 scoring so scoring always reads the latest
cross-organism operon database from depot.
"""
from __future__ import annotations

import argparse
import csv
import fcntl
import sys
from pathlib import Path

_THIS_DIR = Path(__file__).resolve().parent
if str(_THIS_DIR) not in sys.path:
    sys.path.insert(0, str(_THIS_DIR))

import c3_lib as L  # noqa: E402
import c3_occ  # noqa: E402


def _read_organism_from_labeled(labeled_path: Path) -> str:
    with labeled_path.open(newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            return (row.get("organism_name") or "").strip()
    return ""


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", required=True, help="Genome stem / organism name")
    ap.add_argument("--labeled-input", required=True, help="labeling/labeled-genes.tsv")
    ap.add_argument("--operon-info-input", required=True, help="labeling/labeled-genes-operon-info.tsv")
    ap.add_argument("--reference", required=True, help="Depot OCC reference pickle path")
    ap.add_argument("--output-token", required=True, help="Token file marking OCC update completion")
    args = ap.parse_args()

    organism = (args.organism or "").strip()
    labeled_path = Path(args.labeled_input)
    operon_info_path = Path(args.operon_info_input)
    reference_path = Path(args.reference)
    token_path = Path(args.output_token)

    for p in (labeled_path, operon_info_path):
        if not p.is_file():
            print(f"[update-occ-reference] ERROR: input not found: {p}", file=sys.stderr)
            raise SystemExit(1)

    labeled_org = _read_organism_from_labeled(labeled_path)
    if labeled_org and organism and labeled_org != organism:
        print(
            f"[update-occ-reference] ERROR: organism mismatch: --organism={organism} "
            f"but labeled-genes has organism_name={labeled_org}",
            file=sys.stderr,
        )
        raise SystemExit(1)
    if not organism:
        organism = labeled_org
    if not organism:
        print("[update-occ-reference] ERROR: organism could not be determined", file=sys.stderr)
        raise SystemExit(1)

    run_root = labeled_path.parents[2]
    # compute_hash=True -> per-gene aa_hash, so we can fingerprint the genome by
    # CONTENT and record it in the members sidecar (identity for leave-one-out).
    genes = L.load_organism(organism, labeled_path, operon_info_path, compute_hash=True)
    if genes.empty:
        print(f"[update-occ-reference] ERROR: no genes loaded for organism {organism}", file=sys.stderr)
        raise SystemExit(1)

    reference_path.parent.mkdir(parents=True, exist_ok=True)
    token_path.parent.mkdir(parents=True, exist_ok=True)

    lock_path = Path(str(reference_path) + ".lock")
    with lock_path.open("a+") as lock_fh:
        fcntl.flock(lock_fh.fileno(), fcntl.LOCK_EX)
        if reference_path.is_file():
            ref = c3_occ.load_reference(reference_path)
        else:
            ref = c3_occ.new_reference()

        before = len(ref.get("organisms_added", set()))
        already_present = organism in ref.get("organisms_added", set())

        c3_occ.update_reference(ref, genes, run_root, organisms=[organism], skip_existing=True)
        c3_occ.finalize_reference(ref)
        c3_occ.save_reference(ref, reference_path)

        # record this genome's content fingerprint -> token so scoring can find
        # and leave-one-out an already-present genome by content (not by name).
        fingerprint = L.genome_fingerprint(genes["aa_hash"])
        c3_occ.record_member(reference_path, fingerprint, organism)
        # record whole-genome pool stats (gene / operon tallies) for figure
        # provenance -- aggregated (minus the reported genome) at figure time.
        c3_occ.record_pool_stats(reference_path, fingerprint, organism,
                                 L.genome_pool_stats(genes))

        after = len(ref.get("organisms_added", set()))
        added = after - before
        status = "already_present" if already_present else ("added" if added > 0 else "no_change")
        msg = (
            f"occ reference status={status}; organism={organism}; "
            f"fingerprint={fingerprint[:12] or '(none)'}; "
            f"organisms_in_reference={after}; reference={reference_path}"
        )
        token_path.write_text(msg + "\n")

    print(f"[update-occ-reference] {msg}")


if __name__ == "__main__":
    main()
