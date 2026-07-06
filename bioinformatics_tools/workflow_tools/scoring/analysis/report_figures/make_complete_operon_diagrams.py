#!/usr/bin/env python3
"""Complete per-organism operon diagrams.

For ONE organism, draw EVERY multi-gene operon (no named-function gate, nothing
truncated) as a fig06/07/08-style block-arrow map + gene table, grouped by
operon size. Each size bin gets its own file set inside a
`complete-organism-operon-diagrams/` sub-folder of the organism's figures dir:

    complete-organism-operon-diagrams/
        2-gene-operon.tsv            2-gene-operon-p01.png, -p02.png, ...
        3-gene-operon.tsv            3-gene-operon-p01.png, ...
        ...
        79-gene-operon.tsv          79-gene-operon.png

Losslessness: all genes are shown -- operons wider than one row wrap across
stacked rows (arrows keep a constant size, tags run A/B/C… across rows), and the
gene table always lists every member. Pages are split by a height budget so a
size bin with thousands of operons (size 2/3) spans as many pages as needed.

This is an OPT-IN, heavier companion to make_organism_report.py -- run it only
when you want the exhaustive per-organism atlas. Styling is identical to the
report galleries (same reportfig_lib primitives).
"""
from __future__ import annotations
import argparse
import math
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

import reportfig_lib as L

SUBDIR = "complete-organism-operon-diagrams"
_PER_ROW = 10           # genes per row; operons longer than this wrap across rows
                        # (a connecting line links the wrapped rows -- see
                        # draw_gene_track -- so they read as ONE operon)
_TRACK_ROW_IN = L._TRACK_ROW_IN  # figure inches per arrow row (roomy -- wrapped
                        # rows never crowd; arrow thickness is pinned regardless).
                        # Shared with reportfig_lib so the pin scale stays in sync.
_TABLE_FS = 14.0        # gene-table font size (pt). The shared renderer derives the
                        # whole font hierarchy (title = +1, gene tags = +1, etc.)
                        # from this one knob; the line height + page height scale too.
_FIG_W = 22.0           # figure width (in). Wide enough that the bigger-font table
                        # (dynamic descriptor + review-reason columns) reaches the
                        # page margin without crowding (arrows keep a fixed size).
_LINE_IN = round(_TABLE_FS / 72.0 * 1.42, 3)  # inches per gene-table text line,
                        # tied to the font (1.42x line spacing) so rows never
                        # overlap when the font size changes -- ~0.189 in at 9.6pt.
_PAGE_TARGET_IN = 26.0  # soft height budget per page (always >= 1 operon/page)
# Atlas pages are large (up to ~26 in tall) and there are hundreds per organism,
# so full 400-dpi would be ~tens of megapixels each and take many minutes. 170
# dpi keeps the identical fig07/08 styling crisp on screen while cutting raster
# work ~5.5x (cost scales with dpi^2).
_ATLAS_DPI = 170
_DESC_WRAP = 60         # descriptor wrap width (chars): caps the longest descriptor
                        # LINE; the dynamic table layout then sizes the descriptor
                        # column to the actual longest line so names fit on one line.

_RUN_ROOT: Path | None = None
_OPERON_DB: Path | None = None


def _members(op_row) -> list[dict]:
    """All members of an operon, in order -- NEVER capped (lossless)."""
    return L.operon_to_members(op_row)


def _block_height(members: list[dict]) -> float:
    """Figure inches an operon block needs (for pagination): title gap + wrapped
    arrow rows + gap + full table + inter-block gap. Mirrors render_operon_page."""
    n = len(members)
    nrows = math.ceil(n / _PER_ROW) if n > _PER_ROW else 1
    table_in = max(L.member_table_units(members, _DESC_WRAP), 1.0) * _LINE_IN
    return _TRACK_ROW_IN * nrows + table_in + 0.96   # 0.46 title +0.12 intra +0.38 inter


def _paginate(ops: list, heights: list[float]) -> list[list[int]]:
    """Greedily pack operon indices into pages under the height budget; a single
    operon taller than the budget still gets its own page."""
    pages, cur, cur_h = [], [], 0.0
    for i, h in enumerate(heights):
        if cur and cur_h + h > _PAGE_TARGET_IN:
            pages.append(cur)
            cur, cur_h = [], 0.0
        cur.append(i)
        cur_h += h
    if cur:
        pages.append(cur)
    return pages


def _render_page(ops_page: list, size: int, page: int, npages: int,
                 k_lookup, outdir: Path, org_label: str) -> list:
    # Build the operon blocks + their titles, then hand off to the SHARED renderer
    # (reportfig_lib.render_operon_page) so the atlas and the representative-report
    # galleries produce byte-for-byte identical formatting.
    blocks, rows = [], []
    for r in ops_page:
        members = _members(r)
        nrows = math.ceil(len(members) / _PER_ROW) if len(members) > _PER_ROW else 1
        k = k_lookup(r)
        where = f"in {k} pangenome genomes" if k != 1 else "in 1 pangenome genome"
        wrapnote = f"  ({nrows} rows)" if nrows > 1 else ""
        blocks.append({"members": members,
                       "title": f"{size}-gene operon  |  {r['operon_id']}  |  {where}{wrapnote}"})
        rows.append({"operon_id": r["operon_id"], "size": size,
                     "pangenome_organisms": k,
                     "pangenome_occurrences": int(r.get("pangenome_occurrences", 0) or 0),
                     "members_in_order": r["members_in_order"]})
    pagestr = f"  —  page {page}/{npages}" if npages > 1 else ""
    fname = (f"{size}-gene-operon.png" if npages == 1
             else f"{size}-gene-operon-p{page:02d}.png")
    L.render_operon_page(
        outdir / fname, blocks, org_label=org_label,
        suptitle=f"All {size}-gene operons{pagestr}",
        run_root=_RUN_ROOT, note=L.OPERON_CORRECTION_NOTE,
        footer_sources=L.organism_source_lines(_RUN_ROOT, org_label, coords=True,
                                               operon_db=_OPERON_DB),
        fig_width=_FIG_W, table_fs=_TABLE_FS, desc_wrap=_DESC_WRAP,
        per_row=_PER_ROW, min_span=_PER_ROW, dpi=_ATLAS_DPI)
    return rows


def main() -> None:
    global _RUN_ROOT, _OPERON_DB
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-root", required=True)
    ap.add_argument("--organism", required=True)
    ap.add_argument("--operon-db", default=str(L.DEFAULT_OPERON_DB))
    ap.add_argument("--output-dir", default=None,
                    help="defaults to <run>/<organism>/scoring/figures")
    args = ap.parse_args()

    run_root = Path(args.run_root)
    figures = Path(args.output_dir) if args.output_dir else \
        run_root / args.organism / "scoring" / "figures"
    outdir = figures / SUBDIR
    outdir.mkdir(parents=True, exist_ok=True)

    _RUN_ROOT, _OPERON_DB = run_root, Path(args.operon_db)
    L.apply_style()
    org_label = args.organism
    genes = L.load_organism_genes(run_root, args.organism)
    operons = L.build_operons(genes)

    # OCC pool scoping + provenance line, identical to the standard reports:
    # the ACTUAL OCC pool (occ_reference.pkl's organisms), LEAVE-ONE-OUT (this
    # organism was scored against the OTHERS), with gene/operon tallies from the
    # OCC genome-stats sidecar -- NOT discover_organisms(run)+run-folder, which
    # under-counts to just this one genome mid-run (the "1 genome" bug).
    restrict = L.load_occ_organisms() or set(L.discover_organisms(run_root))
    restrict = set(restrict)
    restrict.discard(org_label)
    recurrence = L.load_operon_recurrence(Path(args.operon_db), restrict_to=restrict)
    pool_list = sorted(restrict)
    pstats = L.aggregate_pool_stats(L.load_pool_stats(), pool_list)
    L.set_provenance(L.provenance_text(len(pool_list), pstats, leave_one_out=True))

    ops = operons.copy()
    ops["pangenome_organisms"] = ops["members_in_order"].map(
        lambda m: recurrence.get(m, {}).get("organism_count", 0))
    ops["pangenome_occurrences"] = ops["members_in_order"].map(
        lambda m: recurrence.get(m, {}).get("label_frequency", 0))

    def k_lookup(r):
        return int(r.get("pangenome_organisms", 0) or 0)

    n_pages = 0
    n_operons = 0
    # LARGEST operons first: the big, information-rich operons are the interesting
    # ones, and generating them first stops them being buried behind dozens of
    # pages of 2-/3-gene operons (Haloferax alone: 890 two-gene, ~37 pages).
    for size in sorted(ops["size"].unique(), reverse=True):
        if int(size) < 2:
            continue
        bin_ops = [r for _, r in ops[ops["size"] == size].iterrows()]
        # stable, meaningful order: most-shared operons first, then id
        bin_ops.sort(key=lambda r: (-k_lookup(r), str(r["operon_id"])))
        heights = [_block_height(_members(r)) for r in bin_ops]
        pages = _paginate(bin_ops, heights)
        tsv_rows = []
        for p, idxs in enumerate(pages, start=1):
            page_ops = [bin_ops[i] for i in idxs]
            tsv_rows += _render_page(page_ops, int(size), p, len(pages),
                                     k_lookup, outdir, org_label)
            n_pages += 1
        L.write_tsv(pd.DataFrame(tsv_rows), outdir / f"{int(size)}-gene-operon.tsv")
        n_operons += len(bin_ops)

    print(f"[complete-operon-diagrams] {org_label}: {n_operons} operons across "
          f"{n_pages} pages → {outdir}")


if __name__ == "__main__":
    main()
