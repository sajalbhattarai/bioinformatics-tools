#!/usr/bin/env python3
"""Complete per-organism operon diagrams.

Draws every multi-gene operon of one organism as a block-arrow map with a gene
table (matplotlib, via reportfig_lib), one file set per operon size in
`complete-organism-operon-diagrams/`:

    complete-organism-operon-diagrams/
        2-gene-operon.tsv            2-gene-operon-p01.png, -p02.png, ...
        3-gene-operon.tsv            3-gene-operon-p01.png, ...
        ...
        79-gene-operon.tsv          79-gene-operon.png

Long operons wrap across rows and pages follow a height budget, so every gene
is shown. Optional companion to make_organism_report.py, with the same styling.
"""
from __future__ import annotations
import argparse
import math
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

import reportfig_lib as L

SUBDIR = "complete-organism-operon-diagrams"
_PER_ROW = 10           # genes per row; longer operons wrap across rows
_TRACK_ROW_IN = L._TRACK_ROW_IN  # figure inches per arrow row, shared with reportfig_lib
_TABLE_FS = 14.0        # gene-table font size (pt); the renderer derives other sizes from it
_FIG_W = 22.0           # figure width (in)
_LINE_IN = round(_TABLE_FS / 72.0 * 1.42, 3)  # inches per table line (1.42x line spacing)
_PAGE_TARGET_IN = 26.0  # soft height budget per page (always >= 1 operon/page)
# Lower than the report dpi because atlas pages are large and numerous.
_ATLAS_DPI = 170
_DESC_WRAP = 60         # descriptor wrap width (chars)

_RUN_ROOT: Path | None = None
_OPERON_DB: Path | None = None


def _members(op_row) -> list[dict]:
    """Returns all members of an operon in order."""
    return L.operon_to_members(op_row)


def _block_height(members: list[dict]) -> float:
    """Returns the height in inches of an operon block as render_operon_page draws it (96 px per inch)."""
    return L.operon_block_px(len(members)) / 96.0


def _paginate(ops: list, heights: list[float]) -> list[list[int]]:
    """Packs operon indices greedily into pages under the height budget; an oversize operon gets its own page."""
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
    """Renders one page of operon blocks with reportfig_lib.render_operon_page and returns its TSV rows."""
    blocks, rows = [], []
    for r in ops_page:
        members = _members(r)
        k = k_lookup(r)
        where = f"in {k} pangenome genomes" if k != 1 else "in 1 pangenome genome"
        blocks.append({"members": members, "heading": r["operon_id"],
                       "detail": f"{size}-gene operon  |  {where}"})
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
    """Loads the organism's operons and pangenome recurrence and writes every size bin's pages and TSV."""
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

    # Pool = the OCC reference's organisms minus this one (leave-one-out), with
    # tallies from the OCC genome-stats sidecar, as in the standard reports.
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
        """Returns the number of pangenome organisms sharing the operon."""
        return int(r.get("pangenome_organisms", 0) or 0)

    n_pages = 0
    n_operons = 0
    # Largest operons first, so they are not behind many pages of small ones.
    for size in sorted(ops["size"].unique(), reverse=True):
        if int(size) < 2:
            continue
        bin_ops = [r for _, r in ops[ops["size"] == size].iterrows()]
        # Most-shared operons first, then by id
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
