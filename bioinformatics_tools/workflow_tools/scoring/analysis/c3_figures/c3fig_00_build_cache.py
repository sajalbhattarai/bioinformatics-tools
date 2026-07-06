#!/usr/bin/env python3
"""c3fig_00_build_cache.py — parse all organisms ONCE into pickle caches.

Builds two caches under <output-dir>/_cache/:
  * genes.pkl   — one row per protein-coding gene (all organisms), joined from
                  labeled-genes.tsv + labeled-genes-operon-info.tsv.
  * operons.pkl — one row per operon (derived), with composition, geometry,
                  strand and probability summaries.

Every figure script reads these caches (<1 s) instead of re-parsing ~2 GB TSV.
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import c3fig_lib as L


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--labeling-root", required=True,
                    help="output run dir containing <organism>/labeling/…")
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--no-hash", action="store_true",
                    help="skip AA-sequence hashing (faster)")
    args = ap.parse_args()

    root = Path(args.labeling_root)
    out = Path(args.output_dir)
    cache = out / "_cache"
    cache.mkdir(parents=True, exist_ok=True)

    orgs = L.discover_organisms(root)
    print(f"[build_cache] {len(orgs)} organisms found", file=sys.stderr)

    genes = L.build_gene_table(root, compute_hash=not args.no_hash)
    L.save_cache(genes, cache / "genes.pkl")
    print(f"[build_cache] genes: {len(genes):,} rows", file=sys.stderr)

    operons = L.build_operon_table(genes)
    L.save_cache(operons, cache / "operons.pkl")
    print(f"[build_cache] operons: {len(operons):,} rows", file=sys.stderr)

    pairs = L.build_adjacent_pairs(genes)
    L.save_cache(pairs, cache / "adjacent_pairs.pkl")
    print(f"[build_cache] adjacent_pairs: {len(pairs):,} rows", file=sys.stderr)

    # quick human-readable manifest
    with open(cache / "cache-manifest.txt", "w") as fh:
        fh.write(f"organisms\t{genes['organism'].nunique()}\n")
        fh.write(f"genes_total\t{len(genes)}\n")
        fh.write(f"genes_in_operon\t{int(genes['in_operon'].sum())}\n")
        fh.write(f"genes_uninformative\t{int(genes['uninformative'].sum())}\n")
        fh.write(f"operons_total\t{len(operons)}\n")
        fh.write(f"adjacent_pairs_total\t{len(pairs)}\n")
    print("[build_cache] done", file=sys.stderr)


if __name__ == "__main__":
    main()
