#!/usr/bin/env python3
"""make_circular_genome.py -- a confidence-coloured circular genome map built
ENTIRELY from a FINAL_ANNOTATION_WITH_CONFIDENCE.tsv.

Everything is read from the FINAL table, nothing recomputed except GC% (from the
per-gene nucleotide sequence the table already carries):
  * replicon / contig identity   <- parsed from gene_id  (<accession>_<start><strand><len>)
  * gene arc                      <- RAST_start / RAST_end / RAST_strand
  * ring colour                   <- CONFIDENCE_TIER
  * review track                  <- NEEDS_REVIEW?
  * GC ring                       <- RAST_na_sequence (windowed)

Layout: each replicon (or, for a draft assembly, each contig) is an arc whose
angular width is proportional to its length (with a small floor so tiny
replicons stay visible); a fixed gap wedge separates adjacent arcs. A single
complete chromosome is one near-full circle with one gap at the origin. An
incomplete genome is the SAME picture with more arcs and more gaps -- the gaps
ARE the unassembled/unknown regions.

Usage: make_circular_genome.py <FINAL.tsv> <out.png>
"""
import csv
import glob
import re
import sys
from math import cos, pi, sin
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, PatchCollection
from matplotlib.patches import Wedge

plt.rcParams["font.family"] = "sans-serif"   # Arial/Calibri-style (DejaVu/Liberation Sans)
plt.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "DejaVu Sans"]
csv.field_size_limit(10 ** 8)

# confidence tier -> colour. Okabe-Ito (colour-vision-deficiency safe): the five
# tiers are DISTINCT hues, not a green->red ramp (which collapses under
# deuteranopia). Validated worst-all-pairs CVD deltaE = 16 (target >= 12).
# Confidence tiers are ORDERED (highest > high > medium > fair > low), so they
# get a sequential single-hue ramp, not a categorical rainbow. The previous
# green/blue/yellow/orange/red set failed every ordinal check: lightness was
# non-monotone (L 0.678, 0.679, 0.877, 0.765, 0.670 -- it went up then back
# down), highest and high differed in lightness by 0.001 so they were identical
# in grayscale, the yellow sat at 1.41:1 on white, and the hue spread was 141
# degrees. Rank encoded as hue reads as five unrelated categories.
#
# This ramp is validated (dataviz validate_palette.js --ordinal, light mode):
# monotone lightness, every adjacent gap >= 0.06, 10 degree hue spread, light
# end 2.66:1 against the surface. Darker = higher confidence, so a well
# annotated genome reads solid and a poor one washes out.
#
# Keep this in sync with TIER_COL in gen_genome_viewer.py.
TIER = [
    ("highest", "#0b2842"),   # darkest
    ("high", "#154064"),
    ("medium", "#256291"),
    ("fair", "#4184b5"),
    ("low", "#6ba3c8"),       # lightest
    ("NOT_APPLICABLE_NON_CODING", "#c8c8c8"),
]
TCOL = dict(TIER)
# Review flags are a RESERVED STATUS colour, never a tier step -- an alarm must
# not be confusable with a ranking. Grey worked against the old rainbow only by
# being the one unsaturated thing on the figure; against a single-hue blue ramp
# it reads as just another neutral and the flags stop announcing themselves.
# This red shares no hue with the ramp, so flagged regions are unmistakable.
# Matches FLAG in gen_genome_viewer.py.
REVIEW_COL = "#b32b1e"
GAP_DEG = 4.0                 # angular gap between replicons/contigs (the "unknown" wedge)
GAP_BUDGET_DEG = 40.0         # TOTAL angle spent on gaps, however many contigs there are.
                              # A fixed per-contig gap only works for a finished genome: a
                              # 420-contig draft assembly would want 420*4 = 1680 deg of gap
                              # alone, the sequence budget goes negative, every arc falls to
                              # the floor, and the drawing wraps the circle nine times over.
MIN_SPAN_DEG = 4.0            # floor so a tiny plasmid/contig is still visible
MAX_LABELLED = 12             # past this the map is a SHAPE, not an index -- per-contig
                              # captions become an unreadable wall and the clickable
                              # FINAL_GENOME_VIEWER.html is where identities belong
TICK_MIN_SPAN_DEG = 2.0       # don't hang Mb ticks off a hairline arc
START_DEG = 90.0             # 12 o'clock origin, genome runs clockwise
INK = "#000000"

# ring radii (outer -> inner)
R_BB0, R_BB1 = 1.00, 1.035           # backbone band
R_FWD0, R_FWD1 = 0.905, 0.985        # forward-strand genes
R_REV0, R_REV1 = 0.815, 0.895        # reverse-strand genes
R_REV_TICK0, R_REV_TICK1 = 0.775, 0.805   # review-flag ticks
R_GC = 0.60                          # GC ring baseline
GC_AMP = 0.13                        # GC deviation amplitude
GC_SCALE = 0.09                      # GC deviation (frac) mapped to full amplitude


def col(row, name):
    for k in row:
        if re.sub(r"^Column-[A-Z]+:\s*", "", k or "").strip().lower() == name.lower():
            return row[k] or ""
    return ""


def contig_of(gene_id):
    return re.sub(r"_[0-9]+[+-][0-9]+$", "", gene_id)


def gc(seq):
    seq = seq.upper()
    n = len(seq)
    return (seq.count("G") + seq.count("C")) / n if n else None


def xy(r, deg):
    a = deg * pi / 180.0
    return r * cos(a), r * sin(a)


def main():
    final, out = sys.argv[1], sys.argv[2]
    rows = [r for r in csv.DictReader(open(final, newline=""), delimiter="\t")]
    organism = col(rows[0], "organism_name")

    # group genes by contig, record contig length = max end
    contigs = {}
    for r in rows:
        gid = col(r, "gene_id")
        if not gid:
            continue
        c = contig_of(gid)
        try:
            s, e = int(col(r, "RAST_start")), int(col(r, "RAST_end"))
        except ValueError:
            continue
        d = contigs.setdefault(c, {"len": 0, "genes": []})
        d["len"] = max(d["len"], e)
        d["genes"].append(r)
    order = sorted(contigs, key=lambda c: -contigs[c]["len"])   # largest first
    total_len = sum(contigs[c]["len"] for c in order)
    n_gap = len(order)
    # Gaps share a FIXED total budget, so the sequence budget stays positive at any
    # contig count: a finished genome keeps the full 4 deg wedge, a fragmented draft
    # gets proportionally thinner ones and the gaps still read as ~11% of the circle.
    gap_deg = min(GAP_DEG, GAP_BUDGET_DEG / n_gap)
    span_total = 360.0 - n_gap * gap_deg                        # angle available for sequence

    # angular span per contig, proportional to length but with a floor for tiny ones.
    # The floor is capped at half the per-contig share, so even if EVERY contig is
    # floored the leftover stays positive and the arcs still close at exactly 360.
    min_span = min(MIN_SPAN_DEG, 0.5 * span_total / n_gap)
    floored = [c for c in order if span_total * contigs[c]["len"] / total_len < min_span]
    fixed = len(floored) * min_span
    big_len = sum(contigs[c]["len"] for c in order if c not in floored) or 1
    rest = max(span_total - fixed, 0.0)
    spans = {c: (min_span if c in floored else rest * contigs[c]["len"] / big_len) for c in order}

    # assign each contig a start angle going clockwise from START_DEG
    layout, cur = {}, START_DEG
    for c in order:
        cur -= gap_deg                                          # gap precedes each contig
        layout[c] = (cur, spans[c])                            # start angle (high), span (deg)
        cur -= spans[c]

    def ang(c, pos):
        s0, span = layout[c]
        return s0 - span * (pos / contigs[c]["len"])           # clockwise: decreasing angle

    fig, ax = plt.subplots(figsize=(11, 11))
    fig.patch.set_facecolor("white"); ax.set_facecolor("white")
    ax.set_xlim(-1.36, 1.36); ax.set_ylim(-1.42, 1.30); ax.set_aspect("equal"); ax.axis("off")

    # genome-wide mean GC (length-weighted) for the deviation baseline
    tot_gc = tot_bp = 0.0
    for r in rows:
        seq = col(r, "RAST_na_sequence")
        g = gc(seq)
        if g is not None:
            tot_gc += g * len(seq); tot_bp += len(seq)
    gc_mean = tot_gc / tot_bp if tot_bp else 0.5

    gene_patches, gene_colors, review_lines = [], [], []
    for c in order:
        for r in contigs[c]["genes"]:
            try:
                s, e = int(col(r, "RAST_start")), int(col(r, "RAST_end"))
            except ValueError:
                continue
            a_hi, a_lo = ang(c, s), ang(c, e)
            if a_hi < a_lo:
                a_hi, a_lo = a_lo, a_hi
            tier = col(r, "CONFIDENCE_TIER").strip() or "NOT_APPLICABLE_NON_CODING"
            color = TCOL.get(tier, "#b3b8bf")
            fwd = col(r, "RAST_strand").strip() != "-"
            r0, r1 = (R_FWD0, R_FWD1) if fwd else (R_REV0, R_REV1)
            gene_patches.append(Wedge((0, 0), r1, a_lo, max(a_hi, a_lo + 0.04), width=r1 - r0))
            gene_colors.append(color)
            if col(r, "NEEDS_REVIEW?").strip().lower() == "yes":
                mid = (a_hi + a_lo) / 2
                review_lines.append([xy(R_REV_TICK0, mid), xy(R_REV_TICK1, mid)])

    ax.add_collection(PatchCollection(gene_patches, facecolor=gene_colors, edgecolor="none", antialiased=True))
    if review_lines:
        ax.add_collection(LineCollection(review_lines, colors=REVIEW_COL, linewidths=0.7, alpha=0.9))

    # GC ring: length-weighted mean GC in fixed-bp windows per contig -> smooth filled band
    for c in order:
        L = contigs[c]["len"]
        nb = max(12, min(360, L // 6000))
        wbp = L / nb
        acc = [[0.0, 0.0] for _ in range(nb + 1)]              # [sum gc*bp, sum bp] per window
        for r in contigs[c]["genes"]:
            try:
                s, e = int(col(r, "RAST_start")), int(col(r, "RAST_end"))
            except ValueError:
                continue
            g = gc(col(r, "RAST_na_sequence"))
            if g is None:
                continue
            wi = min(nb, int(((s + e) / 2) / wbp))
            acc[wi][0] += g * (e - s + 1); acc[wi][1] += (e - s + 1)
        curve = [(i, acc[i][0] / acc[i][1]) for i in range(nb + 1) if acc[i][1] > 0]
        if len(curve) < 3:
            continue
        cx, cy = [], []
        for i, m in curve:
            a = ang(c, (i + 0.5) * wbp)
            rr = R_GC + GC_AMP * max(-1.0, min(1.0, (m - gc_mean) / GC_SCALE))
            p = xy(rr, a); cx.append(p[0]); cy.append(p[1])
        bx, by = [], []
        for i, _ in reversed(curve):
            p = xy(R_GC, ang(c, (i + 0.5) * wbp)); bx.append(p[0]); by.append(p[1])
        ax.plot(cx, cy, color="#000000", lw=0.7)          # clean GC line, no shaded fill
    ax.add_patch(plt.Circle((0, 0), R_GC, fill=False, ec="#d5d5d5", lw=0.7))

    # backbone arcs + Mb ticks
    bb = [Wedge((0, 0), R_BB1, layout[c][0] - layout[c][1], layout[c][0], width=R_BB1 - R_BB0) for c in order]
    ax.add_collection(PatchCollection(bb, facecolor="#000000", edgecolor="none"))
    for c in order:
        L = contigs[c]["len"]
        if spans[c] < TICK_MIN_SPAN_DEG:
            continue          # hairline arc: a tick is noise, not a coordinate
        step = 1_000_000 if L > 3_000_000 else 500_000
        # Same threshold as the replicon captions: on a draft assembly every
        # contig restarts its coordinates at zero, so the numerals degenerate
        # into a ring of "0.0" that says nothing. Ticks alone carry the scale
        # there; the centre block carries the total.
        show_nums = L >= 250_000 and len(order) <= MAX_LABELLED
        p = 0
        while p <= L:
            a = ang(c, p)
            x0, y0 = xy(R_BB1, a); x1, y1 = xy(R_BB1 + 0.02, a)
            ax.plot([x0, x1], [y0, y1], color="#000000", lw=0.7)
            if show_nums and p % step == 0:
                xt, yt = xy(R_BB1 + 0.05, a)
                ax.text(xt, yt, f"{p/1e6:.1f}", ha="center", va="center", fontsize=7.5, color="#000000")
            p += step

    # Replicon labels: only for a finished-ish genome (>1 replicon but few enough to
    # caption). A draft assembly's contig names carry no information a reader can use
    # at this scale -- they are looked up by CLICKING the accompanying
    # FINAL_GENOME_VIEWER.html, which carries per-gene and per-operon identity.
    # Printing all of them here buries the map: the anti-collision loop pushes captions
    # out to absurd radii and bbox_inches="tight" then shrinks the figure to a dot.
    if 1 < len(order) <= MAX_LABELLED:
        placed = []                                            # (angle, radius) already used
        for c in order:
            mid = layout[c][0] - layout[c][1] / 2
            rad = R_BB1 + 0.135                                 # clear of the Mb-tick numerals
            while any(abs(((mid - pa + 180) % 360) - 180) < 13 and abs(rad - pr) < 0.075
                      for pa, pr in placed):
                rad += 0.085
            placed.append((mid, rad))
            fs = 8.5 if contigs[c]["len"] >= 250_000 else 7.3
            xl, yl = xy(rad, mid)
            ax.text(xl, yl, f"{c}\n{contigs[c]['len']/1e6:.2f} Mb", ha="center", va="center",
                    fontsize=fs, color="#000000", linespacing=0.95)

    # centre block: stats only (the genome identifier / filename is the figure title)
    n_flag = sum(1 for r in rows if col(r, "NEEDS_REVIEW?").strip().lower() == "yes")
    ax.text(0, 0.05, f"{total_len/1e6:.2f} Mb · {len(rows):,} genes", ha="center", va="center",
            fontsize=13, color="#000000")
    # a genome carved into hundreds of pieces is a draft assembly -- call those contigs,
    # not replicons, since only a handful of them are real replicating molecules
    unit = "contig" if len(order) > MAX_LABELLED else "replicon"
    ax.text(0, -0.05, f"{len(order)} {unit}{'s' if len(order) > 1 else ''} · "
            f"{n_flag:,} flagged for review", ha="center", va="center", fontsize=11.5, color="#000000")

    # tier legend (bottom), only tiers actually present
    used = set(gene_colors)
    present = [(("non-coding" if t.startswith("NOT_") else t), c) for t, c in TIER if c in used]
    widths = [0.05 + 0.022 * len(lab) + 0.055 for lab, _ in present] + [0.05 + 0.022 * len("review flag") + 0.03]
    lx = -sum(widths) / 2
    for lab, c in present:
        ax.add_patch(plt.Rectangle((lx, -1.345), 0.045, 0.045, facecolor=c, edgecolor="none"))
        ax.text(lx + 0.056, -1.322, lab, ha="left", va="center", fontsize=10.5, color="#000000")
        lx += 0.05 + 0.022 * len(lab) + 0.055
    ax.plot([lx + 0.022], [-1.322], marker="|", color=REVIEW_COL, ms=12, mew=1.6)
    ax.text(lx + 0.05, -1.322, "review flag", ha="left", va="center", fontsize=10.5, color="#000000")

    ax.set_title(organism, fontsize=14, color=INK, pad=6)   # genome identifier / filename as the title
    out = Path(out); out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"wrote {out.name}  ({organism}: {len(rows)} genes, {len(order)} {unit}s, {n_flag} flagged)")


if __name__ == "__main__":
    main()
