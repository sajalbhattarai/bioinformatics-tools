#!/usr/bin/env python3
"""Figure 56 - Operon boundaries & partner switching (with a genomic network).

THE QUESTION (from the operon-index thread)
--------------------------------------------
For a candidate gene, are its closest members separated into ANOTHER operon?
If yes, what are the statistics - and when a neighbour IS separated, is it now
partnered with OTHER genes (repartnered) or left alone (a singleton)?  And show
it as a graph.

DEFINITIONS
-----------
* Gene identity = clean_descriptor (functional name).  The SAME descriptor in
  different genomes is the SAME gene.
* Genomic adjacency = two genes consecutive (sorted by start) on the SAME contig,
  regardless of operon membership.  Intergenic distance = down.start-up.end-1.
* An operonic gene = in_operon and member_count >= 2.
* Every genomic adjacency is one of:
      same_operon         both operonic, SAME operon        (a kept partnership)
      diff_operon         both operonic, DIFFERENT operons  (a broken boundary)
      operon_vs_singleton one operonic, one singleton
      both_singleton      neither operonic
* For an operonic gene, a genomic neighbour that is NOT in its operon is
  "separated".  A separated neighbour that is itself operonic has been
  "repartnered" (it is in another operon with other genes); otherwise it is a
  "singleton".

Panels
------
(a) OPERON-BOUNDARY CURVE  P(consecutive genes share an operon | intergenic bin):
    a sharp cliff - below ~100 bp neighbours are usually the same operon, above
    ~100-200 bp they almost never are.  Operon boundaries live in the gap.
(b) classification of every genomic adjacency (counts, %, median gap): kept
    partnerships are tight (~2 bp); boundaries and singleton joins are far apart.
(c) per operonic gene: is your CLOSEST neighbour kept in-operon, and are you an
    interior or a boundary gene?  Closest neighbours are almost always kept.
(d) fate of SEPARATED neighbours: repartnered into another operon vs left a
    singleton.  This answers "are the separated partners now with other genes?"
(e) a real genomic neighbourhood (E. coli) drawn as a network: nodes = genes
    coloured by operon, thick solid links = kept operon partnerships, dashed
    links = broken boundaries; the boundary gap and strand flips are annotated.

Read-only analysis; does NOT modify the scoring pipeline.
"""
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import c3fig_lib as L
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch

csv.field_size_limit(10_000_000)

SUB = "01-operon-context-confidence"

BINS = [(-10**9, -1, "<0"), (0, 0, "0"), (1, 100, "1-100"),
        (101, 200, "100-200"), (201, 300, "200-300"), (301, 400, "300-400"),
        (401, 500, "400-500"), (501, 10**9, ">500")]
BIN_NAMES = [b[2] for b in BINS]
BIN_COLORS = [L.GREEN, L.TEAL, L.LIME, L.YELLOW, L.AMBER, L.ORANGE, L.RED, "#8a1220"]

CLASS_ORDER = ["same_operon", "diff_operon", "operon_vs_singleton", "both_singleton"]
CLASS_COLOR = {"same_operon": L.BLUE, "diff_operon": L.RED,
               "operon_vs_singleton": L.ORANGE, "both_singleton": "#9e9e9e"}
CLASS_NAME = {"same_operon": "same operon\n(kept partner)",
              "diff_operon": "different operons\n(broken boundary)",
              "operon_vs_singleton": "operon +\nsingleton",
              "both_singleton": "both\nsingletons"}


def _bin_of(gap):
    for lo, hi, nm in BINS:
        if lo <= gap <= hi:
            return nm
    return ">500"


def _build_contig_map(run_root, organisms):
    rx = re.compile(r"^(.*)_(\d+)([+-])(\d+)$")
    m = {}
    for org in organisms:
        path = run_root / org / "labeling" / "labeled-genes.tsv"
        if not path.is_file():
            continue
        with open(path, newline="") as fh:
            rd = csv.reader(fh, delimiter="\t")
            header = next(rd)
            iF = header.index("feature_id")
            iG = header.index("gene_id")
            for row in rd:
                if len(row) > max(iF, iG):
                    hit = rx.match(row[iG] or "")
                    m[(org, row[iF])] = hit.group(1) if hit else None
    return m


def make(genes, operons, outdir):
    run_root = outdir.parents[3]
    g = genes.copy()
    organisms = sorted(g["organism"].unique())
    contig = _build_contig_map(run_root, organisms)
    g["contig"] = [contig.get((o, f)) for o, f in zip(g["organism"], g["feature_id"])]

    inop = {(o, f): bool(io) for o, f, io in
            zip(g["organism"], g["feature_id"], g["in_operon"])}
    mcnt = {(o, f): int(m) for o, f, m in
            zip(g["organism"], g["feature_id"], g["member_count"])}

    def is_op(o, f):
        return inop.get((o, f), False) and mcnt.get((o, f), 0) >= 2

    # ---- walk every contig: genomic adjacencies + per-gene boundary status --
    cls_count = defaultdict(int)
    cls_gap = defaultdict(list)
    bin_same = defaultdict(lambda: [0, 0])   # bin -> [n_same_operon, n_total]
    sep_repartnered = 0
    sep_singleton = 0
    closest_inop = 0
    closest_sep = 0
    interior = 0
    boundary = 0
    total_op = 0

    for (org, ctg), sub in g.groupby(["organism", "contig"]):
        sub = sub.sort_values("start")
        r = sub.to_dict("records")
        for a, b in zip(r[:-1], r[1:]):
            gap = int(b["start"]) - int(a["end"]) - 1
            ao, bo = is_op(org, a["feature_id"]), is_op(org, b["feature_id"])
            if ao and bo and a["operon_id"] == b["operon_id"]:
                cls = "same_operon"
            elif ao and bo:
                cls = "diff_operon"
            elif ao or bo:
                cls = "operon_vs_singleton"
            else:
                cls = "both_singleton"
            cls_count[cls] += 1
            cls_gap[cls].append(gap)
            nm = _bin_of(gap)
            bin_same[nm][1] += 1
            if cls == "same_operon":
                bin_same[nm][0] += 1
        # per operonic gene: closest / interior-boundary / fate of separated
        for i, gene in enumerate(r):
            if not is_op(org, gene["feature_id"]):
                continue
            total_op += 1
            nbrs = []
            if i > 0:
                nbrs.append((r[i - 1], int(gene["start"]) - int(r[i - 1]["end"]) - 1))
            if i < len(r) - 1:
                nbrs.append((r[i + 1], int(r[i + 1]["start"]) - int(gene["end"]) - 1))
            sep_here = False
            for nb, gap in nbrs:
                same = is_op(org, nb["feature_id"]) and nb["operon_id"] == gene["operon_id"]
                if same:
                    continue
                sep_here = True
                if is_op(org, nb["feature_id"]):
                    sep_repartnered += 1
                else:
                    sep_singleton += 1
            if nbrs:
                nb, gap = min(nbrs, key=lambda t: t[1])
                same = is_op(org, nb["feature_id"]) and nb["operon_id"] == gene["operon_id"]
                if same:
                    closest_inop += 1
                else:
                    closest_sep += 1
            if sep_here:
                boundary += 1
            else:
                interior += 1

    n_adj = sum(cls_count.values()) or 1

    # ======================= FIGURE =========================================
    fig = plt.figure(figsize=(15.8, 16.2))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 1.0, 1.08],
                          hspace=0.42, wspace=0.24)
    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1])
    axC = fig.add_subplot(gs[1, 0])
    axD = fig.add_subplot(gs[1, 1])
    axE = fig.add_subplot(gs[2, :])

    # ---- (a) operon-boundary curve ----------------------------------------
    xs, ps, ns = [], [], []
    for nm in BIN_NAMES:
        same, tot = bin_same.get(nm, [0, 0])
        if tot:
            xs.append(nm)
            ps.append(100 * same / tot)
            ns.append(tot)
    x = np.arange(len(xs))
    cols = [BIN_COLORS[BIN_NAMES.index(nm)] for nm in xs]
    axA.bar(x, ps, color=cols, edgecolor="black", linewidth=0.6, width=0.74)
    for xi, p, nn in zip(x, ps, ns):
        axA.annotate("%.1f%%\nn=%d" % (p, nn), (xi, p), ha="center", va="bottom",
                     fontsize=8.2, fontweight="bold")
    axA.set_xticks(x)
    axA.set_xticklabels(xs, rotation=25, ha="right")
    axA.set_ylim(0, 100)
    axA.set_xlabel("Intergenic distance between genomic neighbours (bp)",
                   fontweight="bold")
    axA.set_ylabel("P(neighbours share an operon)  %", fontweight="bold")
    axA.set_title("(a) The operon-boundary cliff\n"
                  "<=100 bp -> usually one operon;  >100 bp -> almost never")
    L.boldticks(axA)
    axA.grid(False)

    # ---- (b) classification of genomic adjacencies ------------------------
    vals = [cls_count.get(c, 0) for c in CLASS_ORDER]
    cols = [CLASS_COLOR[c] for c in CLASS_ORDER]
    names = [CLASS_NAME[c] for c in CLASS_ORDER]
    xb = np.arange(len(CLASS_ORDER))
    axB.bar(xb, vals, color=cols, edgecolor="black", linewidth=0.6, width=0.72)
    for xi, c in zip(xb, CLASS_ORDER):
        v = cls_count.get(c, 0)
        med = int(np.median(cls_gap[c])) if cls_gap[c] else 0
        axB.annotate("%d\n%.1f%%\nmed %dbp" % (v, 100 * v / n_adj, med),
                     (xi, v), ha="center", va="bottom", fontsize=8.0,
                     fontweight="bold")
    axB.set_xticks(xb)
    axB.set_xticklabels(names, fontsize=8.5)
    axB.set_ylim(0, max(vals) * 1.24 if max(vals) else 1)
    axB.set_ylabel("Genomic adjacencies", fontweight="bold")
    axB.set_title("(b) Every neighbour relationship, classified\n"
                  "kept partnerships are tight; boundaries sit in wide gaps")
    L.boldticks(axB)
    axB.grid(False)

    # ---- (c) per-gene: closest neighbour + boundary status ----------------
    c_labels = ["closest nbr\nKEPT in-operon", "closest nbr\nseparated",
                "interior gene\n(both kept)", "boundary gene\n(>=1 separated)"]
    c_vals = [100 * closest_inop / total_op, 100 * closest_sep / total_op,
              100 * interior / total_op, 100 * boundary / total_op]
    c_raw = [closest_inop, closest_sep, interior, boundary]
    c_cols = [L.BLUE, L.RED, L.CYAN, L.ORANGE]
    xc = np.arange(len(c_labels))
    axC.bar(xc, c_vals, color=c_cols, edgecolor="black", linewidth=0.6, width=0.7)
    for xi, p, nn in zip(xc, c_vals, c_raw):
        axC.annotate("%.1f%%\n(%d)" % (p, nn), (xi, p), ha="center", va="bottom",
                     fontsize=8.4, fontweight="bold")
    axC.set_xticks(xc)
    axC.set_xticklabels(c_labels, fontsize=8.5)
    axC.set_ylim(0, 108)
    axC.set_ylabel("Operonic genes  %", fontweight="bold")
    axC.set_title("(c) Are a gene's CLOSEST members separated?\n"
                  "no - %.0f%% keep their nearest neighbour; boundaries are on the far side"
                  % (100 * closest_inop / total_op))
    L.boldticks(axC)
    axC.grid(False)

    # ---- (d) fate of separated neighbours ---------------------------------
    n_sep = sep_repartnered + sep_singleton
    d_vals = [sep_repartnered, sep_singleton]
    d_pct = [100 * sep_repartnered / n_sep, 100 * sep_singleton / n_sep]
    d_names = ["repartnered\n(in ANOTHER operon\nwith other genes)",
               "singleton\n(alone, no operon)"]
    d_cols = [L.PURPLE, "#9e9e9e"]
    xd = np.arange(2)
    axD.bar(xd, d_vals, color=d_cols, edgecolor="black", linewidth=0.6, width=0.62)
    for xi, v, p in zip(xd, d_vals, d_pct):
        axD.annotate("%d\n%.1f%%" % (v, p), (xi, v), ha="center", va="bottom",
                     fontsize=9.5, fontweight="bold")
    axD.set_xticks(xd)
    axD.set_xticklabels(d_names, fontsize=9)
    axD.set_ylim(0, max(d_vals) * 1.2 if max(d_vals) else 1)
    axD.set_ylabel("Separated genomic neighbours", fontweight="bold")
    axD.set_title("(d) When a neighbour IS separated, where does it go?\n"
                  "%d separated: %.0f%% repartnered, %.0f%% left singleton"
                  % (n_sep, d_pct[0], d_pct[1]))
    L.boldticks(axD)
    axD.grid(False)

    # ---- (e) genomic-neighbourhood network --------------------------------
    _draw_network(axE, g, is_op, organisms)

    fig.suptitle("Operon boundaries and partner switching: neighbours are kept "
                 "in-operon across tight gaps, separated across wide ones",
                 fontsize=13.5, fontweight="bold", y=0.997)

    # ---- TSVs --------------------------------------------------------------
    bcurve = pd.DataFrame([{
        "intergenic_bin": nm, "n_adjacencies": bin_same.get(nm, [0, 0])[1],
        "n_same_operon": bin_same.get(nm, [0, 0])[0],
        "p_same_operon": (bin_same.get(nm, [0, 0])[0] / bin_same.get(nm, [0, 0])[1])
        if bin_same.get(nm, [0, 0])[1] else float("nan"),
    } for nm in BIN_NAMES])
    L.write_tsv(bcurve, outdir / "fig02_boundary_curve.tsv")

    clsdf = pd.DataFrame([{
        "class": c, "n_adjacencies": cls_count.get(c, 0),
        "pct": round(100 * cls_count.get(c, 0) / n_adj, 3),
        "median_gap_bp": int(np.median(cls_gap[c])) if cls_gap[c] else 0,
        "mean_gap_bp": round(float(np.mean(cls_gap[c])), 2) if cls_gap[c] else 0.0,
    } for c in CLASS_ORDER])
    L.write_tsv(clsdf, outdir / "fig02_adjacency_classes.tsv")

    summ = pd.DataFrame([{
        "n_genomic_adjacencies": n_adj,
        "n_operonic_genes": total_op,
        "pct_closest_neighbour_kept": round(100 * closest_inop / total_op, 2),
        "pct_closest_neighbour_separated": round(100 * closest_sep / total_op, 2),
        "pct_interior_genes": round(100 * interior / total_op, 2),
        "pct_boundary_genes": round(100 * boundary / total_op, 2),
        "n_separated_neighbours": n_sep,
        "pct_separated_repartnered": round(100 * sep_repartnered / n_sep, 2),
        "pct_separated_singleton": round(100 * sep_singleton / n_sep, 2),
    }])
    L.write_tsv(summ, outdir / "fig02_summary.tsv")

    L.savefig(fig, outdir / "fig02_operon_boundaries_partner_switching.png")


def _find_window(g, is_op, organisms):
    """First clean E. coli boundary (operon>=3 | gap>=150 | operon>=2), expanded
    to whole operons, capped ~11 genes.  Returns list of gene records."""
    ec = [o for o in organisms if o.startswith("Escherichia")]
    ec = ec[0] if ec else organisms[0]
    sub = g[g["organism"] == ec]
    ctg = sub["contig"].value_counts().index[0]
    recs = sub[sub["contig"] == ctg].sort_values("start").to_dict("records")

    def opid(rec):
        oid = rec["operon_id"]
        return oid if str(oid).startswith("operon_") and is_op(ec, rec["feature_id"]) else None

    counts = defaultdict(int)
    for rec in recs:
        oi = opid(rec)
        if oi:
            counts[oi] += 1
    for i in range(len(recs) - 1):
        a, b = recs[i], recs[i + 1]
        oa, ob = opid(a), opid(b)
        gap = int(b["start"]) - int(a["end"]) - 1
        if oa and ob and oa != ob and gap >= 150 and counts[oa] >= 3 and counts[ob] >= 2:
            lo, hi = max(0, i - 4), min(len(recs), i + 6)
            # expand to whole operons at the edges
            while lo > 0 and opid(recs[lo]) and opid(recs[lo]) == opid(recs[lo - 1]):
                lo -= 1
            while hi < len(recs) and opid(recs[hi - 1]) and opid(recs[hi - 1]) == opid(recs[hi] if hi < len(recs) else recs[hi - 1]):
                hi += 1
            return ec, recs[lo:hi]
    return ec, recs[:11]


def _draw_network(ax, g, is_op, organisms):
    ec, win = _find_window(g, is_op, organisms)

    def opid(rec):
        oid = rec["operon_id"]
        return oid if str(oid).startswith("operon_") and is_op(ec, rec["feature_id"]) else None

    op_ids = [o for o in dict.fromkeys(opid(r) for r in win) if o]
    op_col = {o: L.BRIGHT[i % len(L.BRIGHT)] for i, o in enumerate(op_ids)}
    SING = "#c8c8c8"

    n = len(win)
    xs = list(range(n))
    ys = [0.45 if r["strand"] == "+" else -0.45 for r in win]

    # operon grouping bands
    j = 0
    while j < n:
        oi = opid(win[j])
        if oi:
            k = j
            while k + 1 < n and opid(win[k + 1]) == oi:
                k += 1
            x0, x1 = xs[j] - 0.42, xs[k] + 0.42
            band = FancyBboxPatch((x0, -1.02), x1 - x0, 2.04,
                                  boxstyle="round,pad=0.02,rounding_size=0.12",
                                  linewidth=1.4, edgecolor=op_col[oi],
                                  facecolor=op_col[oi], alpha=0.12, zorder=0)
            ax.add_patch(band)
            ax.text((x0 + x1) / 2, 1.16, oi.replace("operon_", "op "),
                    ha="center", va="bottom", fontsize=8.5, fontweight="bold",
                    color=op_col[oi])
            j = k + 1
        else:
            j += 1

    # edges between consecutive genes
    for i in range(n - 1):
        a, b = win[i], win[i + 1]
        gap = int(b["start"]) - int(a["end"]) - 1
        same = opid(a) is not None and opid(a) == opid(b)
        if same:
            ax.plot([xs[i], xs[i + 1]], [ys[i], ys[i + 1]], color=op_col[opid(a)],
                    lw=3.2, solid_capstyle="round", zorder=1)
        else:
            ax.plot([xs[i], xs[i + 1]], [ys[i], ys[i + 1]], color="black",
                    lw=1.3, ls=(0, (4, 3)), zorder=1)
        ax.annotate("%d bp" % gap, ((xs[i] + xs[i + 1]) / 2, (ys[i] + ys[i + 1]) / 2),
                    xytext=(0, 7 if not same else 0), textcoords="offset points",
                    ha="center", va="bottom", fontsize=7.2,
                    fontweight="bold", color="black" if not same else "#333333")

    # nodes
    for i, r in enumerate(win):
        oi = opid(r)
        col = op_col[oi] if oi else SING
        ax.scatter([xs[i]], [ys[i]], s=780, c=col, edgecolors="black",
                   linewidths=1.1, zorder=3)
        arrow = ">" if r["strand"] == "+" else "<"
        ax.text(xs[i], ys[i], arrow, ha="center", va="center", fontsize=11,
                fontweight="bold", color="white", zorder=4)
        ax.text(xs[i], ys[i] + (0.30 if ys[i] > 0 else -0.30),
                L.short_desc(r["clean_descriptor"], 16),
                ha="center", va="bottom" if ys[i] > 0 else "top",
                fontsize=7.4, fontweight="bold", rotation=0)

    from matplotlib.lines import Line2D
    handles = [
        Line2D([0], [0], color="black", lw=3.2, label="kept operon partnership"),
        Line2D([0], [0], color="black", lw=1.3, ls=(0, (4, 3)),
               label="broken boundary (separated)"),
        Line2D([0], [0], marker="o", linestyle="none", markersize=10,
               markerfacecolor=SING, markeredgecolor="black",
               label="singleton (no operon)"),
    ]
    ax.legend(handles=handles, loc="lower center", ncol=3, frameon=True,
              fontsize=8.8, bbox_to_anchor=(0.5, -0.02))
    ax.set_xlim(-0.8, n - 0.2)
    ax.set_ylim(-1.5, 1.5)
    ax.set_yticks([0.45, -0.45])
    ax.set_yticklabels(["+ strand", "- strand"], fontweight="bold")
    ax.set_xticks([])
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.set_title("(e) A real genomic neighbourhood (%s): kept partners sit in tight, "
                 "same-operon runs; wide gaps / strand flips break the boundary and "
                 "the separated gene is repartnered into its own operon"
                 % L.short_label(ec), fontsize=10.5, pad=14)
    ax.grid(False)


if __name__ == "__main__":
    L.figure_main(make, SUB)
