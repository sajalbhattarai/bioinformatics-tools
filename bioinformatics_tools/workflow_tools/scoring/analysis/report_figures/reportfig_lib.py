#!/usr/bin/env python3
"""reportfig_lib.py -- shared foundation for the post-scoring report figures.

INDEPENDENT of the scoring pipeline. It reads only FINISHED scoring outputs
(per organism, under a timestamped run folder) plus the depot pangenome operon
reference. It never writes into any scoring folder and cannot affect scoring.

Design rules (enforced by convention across every figure that uses this lib):
  * PRESENTATION ONLY. Figures report the numbers; titles/labels describe what
    is plotted, never a conclusion or judgement (those go stale as the database
    grows). No "this proves / fully assembles / is real" text anywhere.
  * PLAIN LANGUAGE. No jargon in anything a viewer reads (no "cliff / module /
    flagship / link floor / fragmentation / pool of trustable modules"). Use
    "separation of operon members", "operon-context score", etc.
  * EVERY figure emits a .png (>=400 dpi) AND a companion .tsv holding the exact
    plotted numbers, so the figure is auditable from its TSV alone.
  * Colorblind-safe, vibrant palette (Okabe-Ito); panel letters sit ABOVE the
    plot area, clear of the title; no field/text overlaps.
"""
from __future__ import annotations

import re
import sys
import textwrap
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # headless / SLURM
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyBboxPatch, Polygon  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

# Reuse the validated, pipeline-canonical helpers (read-only) for descriptor
# cleaning and the uninformative gate, so "gene identity" here matches scoring.
_C3_DIR = Path(__file__).resolve().parent.parent  # .../scoring/analysis
_SCORING_DIR = _C3_DIR.parent                      # .../scoring
for _p in (str(_SCORING_DIR), str(_C3_DIR / "c3_figures")):
    if _p not in sys.path:
        sys.path.insert(0, _p)
try:
    import c3_lib as _c3lib  # clean_descriptor, is_uninformative
    clean_descriptor = _c3lib.clean_descriptor
    is_uninformative = _c3lib.is_uninformative
except Exception:  # pragma: no cover - fall back to a minimal cleaner
    def clean_descriptor(desc: str) -> str:
        d = (desc or "").split("##")[0].strip()
        return re.sub(r"^[A-Za-z0-9 _()-]+:\s*", "", d).strip()

    def is_uninformative(desc: str) -> bool:
        d = (desc or "").strip().lower()
        return (not d) or d in {"hypothetical protein", "unknown", "uncharacterized protein"}

DPI = 400

# ---------------------------------------------------------------------------
# Palette -- vibrant, saturated hues (the house palette from the reference C3
# figures). Categorical hues assigned in fixed order, never cycled past it.
# ---------------------------------------------------------------------------
BLUE = "#1f77ff"
ORANGE = "#ff8c00"
GREEN = "#00b84d"
RED = "#ee2233"
PURPLE = "#9b30ff"
CYAN = "#12c4e6"
YELLOW = "#ffcc00"
AMBER = "#ffb200"
TEAL = "#00b3a4"
PINK = "#ff4da6"
LIME = "#8ce65a"
VERMILLION = RED   # aliases kept so existing figure code stays valid
SKY = CYAN
GREY = "#7a7a7a"
LIGHTGREY = "#c2c2c2"
INK = "#111111"
_DETAIL_PURPLE = "#6f42c1"   # informative detail text (provenance, table sub-headers)
_MAGENTA = "#c71585"         # sources footer ("Files used to build this figure")
OLIVE = "#6b8e23"     # "raised by operon context"
DARKRED = "#8b0000"   # "lowered by operon context"
RAISED, LOWERED = OLIVE, DARKRED
CATEGORICAL = [BLUE, ORANGE, GREEN, RED, PURPLE, CYAN, AMBER, TEAL, PINK, LIME]

# Confidence tiers: fixed order, vibrant "traffic-light" ramp (best->worst) so
# a tier is always the same hue everywhere. NON_CODING kept separate in grey.
CONF_TIER_ORDER = ["highest", "high", "medium", "fair", "low"]
CONF_TIER_COLOR = {
    "highest": "#1f77ff",
    "high": "#00b84d",
    "medium": "#ffcc00",
    "fair": "#ff8c00",
    "low": "#ee2233",
}
NONCODING_TIER = "NOT_APPLICABLE_NON_CODING"
NONCODING_COLOR = LIGHTGREY

# Component semantic colors (kept consistent across all figures).
COMPONENT_COLOR = {
    "C1": BLUE, "C2": ORANGE, "C3": GREEN, "C4": PURPLE,
    "preliminary": CYAN, "final": "#1f4fff",
}
COMPONENT_LABEL = {
    "C1": "C1  database coverage",
    "C2": "C2  operon probability",
    "C3": "C3  operon-context score",
    "C4": "C4  EC-number agreement",
}

# Serif preference: real Times New Roman on machines that have it (production),
# metric-compatible substitutes next, DejaVu Serif as the always-present floor.
_SERIF = ["Times New Roman", "Nimbus Roman No9 L", "Nimbus Roman",
          "Liberation Serif", "Tinos", "Times", "DejaVu Serif"]

# Journal-standard sans-serif: Nimbus Sans is URW's Helvetica-metric clone (i.e.
# Helvetica), Liberation Sans is the Arial clone -- both render identically to
# Helvetica/Arial for print. Preference order, with DejaVu Sans as last resort.
_SANS = ["Helvetica", "Nimbus Sans", "Arial", "Liberation Sans",
         "TeX Gyre Heros", "DejaVu Sans"]

# Font files to register with matplotlib if present -- its default cache omits the
# urw-base35 / liberation-sans dirs, so "Helvetica"/"Nimbus Sans" would otherwise
# silently fall back to DejaVu Sans.
_FONT_FILES = [
    "/usr/share/fonts/urw-base35/NimbusSans-Regular.otf",
    "/usr/share/fonts/urw-base35/NimbusSans-Bold.otf",
    "/usr/share/fonts/urw-base35/NimbusSans-Italic.otf",
    "/usr/share/fonts/urw-base35/NimbusSans-BoldItalic.otf",
    "/usr/share/fonts/liberation-sans/LiberationSans-Regular.ttf",
    "/usr/share/fonts/liberation-sans/LiberationSans-Bold.ttf",
    "/usr/share/fonts/liberation-sans/LiberationSans-Italic.ttf",
    # serif (Times-metric): Nimbus Roman = URW's Times clone; Liberation Serif next.
    "/usr/share/fonts/urw-base35/NimbusRoman-Regular.otf",
    "/usr/share/fonts/urw-base35/NimbusRoman-Bold.otf",
    "/usr/share/fonts/urw-base35/NimbusRoman-Italic.otf",
    "/usr/share/fonts/urw-base35/NimbusRoman-BoldItalic.otf",
    "/usr/share/fonts/liberation-serif/LiberationSerif-Regular.ttf",
    "/usr/share/fonts/liberation-serif/LiberationSerif-Bold.ttf",
    "/usr/share/fonts/liberation-serif/LiberationSerif-Italic.ttf",
    "/usr/share/fonts/liberation-serif/LiberationSerif-BoldItalic.ttf",
]

_STYLE = {
    "font.family": "serif",
    "font.sans-serif": _SANS,
    "font.serif": _SERIF,
    "font.weight": "bold",
    "font.size": 13,
    "axes.titlesize": 15,
    "axes.titleweight": "bold",
    "axes.labelsize": 14,
    "axes.labelweight": "bold",
    "axes.edgecolor": "#222222",
    "axes.linewidth": 1.1,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "text.color": INK,
    "axes.labelcolor": INK,
    "axes.titlecolor": INK,
    "legend.fontsize": 12,
    "legend.frameon": False,
    "figure.titlesize": 17,
    "figure.titleweight": "bold",
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.dpi": DPI,
    "savefig.facecolor": "white",
    "axes.grid": False,          # no grids anywhere
    # mathtext used only to italicise scientific names; force a NON-bold SERIF
    # italic so the binomial is italic but never bold, matching the serif body.
    "mathtext.fontset": "custom",
    "mathtext.rm": "serif",
    "mathtext.it": "serif:italic",
    "mathtext.bf": "serif:bold",
    "mathtext.default": "regular",
}


def _register_fonts() -> None:
    """Register the Helvetica-clone font files so font.family actually resolves
    (matplotlib's default cache omits these dirs). Silently no-ops on any missing
    file / error, leaving the DejaVu Sans fallback."""
    from matplotlib import font_manager as fm
    for p in _FONT_FILES:
        try:
            if Path(p).exists():
                fm.fontManager.addfont(p)
        except Exception:
            pass


def apply_style() -> None:
    _register_fonts()
    plt.rcParams.update(_STYLE)


# ---------------------------------------------------------------------------
# Layout helpers -- panel letter sits ABOVE the plot area and BELOW the title,
# at the top-left, never overlapping the title, the plot, or facet headers.
# ---------------------------------------------------------------------------
def panel_letter(ax, letter: str) -> None:
    # sits at the far top-left corner, below the centred panel title; placed
    # well to the left so a wide centred title never reaches down onto it
    ax.text(-0.065, 1.012, f"({letter})", transform=ax.transAxes,
            fontsize=15, fontweight="bold", va="bottom", ha="left",
            color=INK, clip_on=False)


def set_title(ax, text: str) -> None:
    """Bold descriptive title, padded so it sits clearly ABOVE the panel letter
    (which sits just above the plot)."""
    ax.set_title(text, pad=26, loc="center", fontweight="bold")


# ---------------------------------------------------------------------------
# Provenance line -- a small parenthetical under the title recording the OCC /
# operon-database pool this figure was generated against (organism count + total
# genes). Set once per driver run with set_provenance(); finish() and the two
# raw-suptitle figures then draw it automatically so EVERY figure of a given
# generation carries the same traceable pool size.
_PROVENANCE: str | None = None


def set_provenance(text: str | None) -> None:
    """Register the provenance line drawn under every figure title this run."""
    global _PROVENANCE
    _PROVENANCE = text


def provenance_text(pool_organisms: int, stats=None, leave_one_out: bool = False) -> str:
    """Provenance string: the OCC/operon-database pool this generation was scored
    against. ``stats`` is an aggregate_pool_stats() dict (genes + operon tallies);
    for backward compatibility a bare int is treated as total_genes only.
    ``leave_one_out`` appends a note that THIS organism was excluded from the pool
    (per-organism reports), so the counts carry no candidate-self bias."""
    loo = (" | this candidate genome EXCLUDED from the pool (leave-one-out, "
           "avoids candidate-self bias)" if leave_one_out else "")
    if isinstance(stats, dict):
        return (f"(operon-database pool: {pool_organisms:,} genomes | "
                f"{stats.get('total_genes', 0):,} genes | "
                f"{stats.get('n_operons', 0):,} operons "
                f"[{stats.get('n_informative_operons', 0):,} informative / "
                f"{stats.get('n_uninformative_operons', 0):,} uninformative] | "
                f"{stats.get('singleton_genes', 0):,} non-operonic genes{loo})")
    return (f"(operon-database pool: {pool_organisms:,} genomes | "
            f"{(stats or 0):,} genes{loo})")


def draw_provenance_line(fig, y: float, fontsize: float = 9.3):
    """Draw the provenance line at figure-fraction height `y` if one is set.
    Returns the artist (or None) so the caller can add it to the crop bbox."""
    if not _PROVENANCE:
        return None
    return fig.text(0.5, y, _PROVENANCE, ha="center", va="top", fontsize=fontsize,
                    color=_DETAIL_PURPLE, style="italic")


def pool_total_genes(run_root, organisms) -> int:
    """Total gene count across the given organisms' scored output (the OCC pool).
    Reads each organism's final annotation once; missing organisms are skipped."""
    total = 0
    for org in organisms:
        try:
            total += len(load_organism_genes(run_root, org))
        except Exception:
            pass
    return total


def finish(fig, suptitle: str | None = None, organism: str | None = None,
           top: float = 0.93, h_pad: float = 2.6, w_pad: float = 2.8,
           band: float | None = None) -> None:
    """tight_layout that reserves headroom for the title. When `organism` is
    given, a bold description and a SECOND, non-bold italic organism line are
    drawn in a reserved title band of `band` inches -- big enough that panels
    carrying their own titles never reach up into the organism line."""
    if not suptitle and not organism:
        fig.tight_layout(h_pad=h_pad, w_pad=w_pad)
        return
    # tight_layout(rect=) is unreliable when a gene-track axis is present (it
    # warns and ignores the reserved band), so lay panels out first, THEN
    # reserve a fixed inch title band with subplots_adjust (deterministic). A
    # bigger band is used when an organism line is present (two title lines) or
    # by request; otherwise a single suptitle needs less headroom.
    fh = max(fig.get_figheight(), 3.0)
    b = band if band is not None else (1.55 if organism else 1.4)
    if _PROVENANCE:
        b += 0.30  # extra headroom for the provenance line under the title
    try:
        fig.tight_layout(h_pad=h_pad, w_pad=w_pad)
    except Exception:
        pass
    fig.subplots_adjust(top=1 - b / fh)
    st = fig.suptitle(suptitle or "", y=1 - 0.42 / fh, fontweight="bold")
    extras = [st]
    prov_y = 0.82  # inches-from-top for the provenance line (single-title case)
    if organism:
        draw_organism_line(fig, organism, y=1 - 0.82 / fh, fontsize=13.5)
        prov_y = 1.14
    pv = draw_provenance_line(fig, 1 - prov_y / fh)
    if pv is not None:
        extras.append(pv)
    if organism or pv is not None:
        fig._report_extra_artists = getattr(fig, "_report_extra_artists", []) + extras


def savefig(fig, path: Path, dpi: int = DPI) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    extra = getattr(fig, "_report_extra_artists", None)
    if extra:
        # The title lines (suptitle + organism) and the sources footer live in
        # reserved bands outside the axes. bbox_inches="tight" RE-RENDERS and
        # remaps figure-fraction artists onto the axes; instead compute the tight
        # bbox once (including those artists AND the suptitle, which
        # get_tightbbox otherwise drops when bbox_extra_artists is passed) and
        # crop to that explicit window -- no remap.
        r = fig.canvas.get_renderer()
        fig.draw(r)
        arts = list(extra)
        st = getattr(fig, "_suptitle", None)
        if st is not None and st not in arts:
            arts.append(st)
        # pad the tight bbox (an explicit Bbox gets NO margin, unlike
        # bbox_inches="tight" which pads 0.1") so nothing is flush/clipped at the
        # edges when the PNG is downloaded
        bb = fig.get_tightbbox(r, bbox_extra_artists=arts).padded(0.40)
        fig.savefig(path, dpi=dpi, bbox_inches=bb, facecolor="white")
    else:
        fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"[reportfig] wrote {path.name}", file=sys.stderr)


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)
    print(f"[reportfig] wrote {path.name}", file=sys.stderr)


def short_desc(desc: str, maxlen: int = 26) -> str:
    d = clean_descriptor(desc) or (desc or "")
    d = re.sub(r"\s+", " ", d).strip()
    return d if len(d) <= maxlen else d[: maxlen - 1] + "…"


def short_organism(org: str) -> str:
    """'Escherichia_coli_str._K-12...GCF_x' -> 'E. coli str. K-12'."""
    name = re.sub(r"_GCF_.*$", "", org or "").replace("_", " ").strip()
    parts = name.split()
    if len(parts) >= 2 and parts[0][:1].isupper():
        return f"{parts[0][0]}. " + " ".join(parts[1:4])
    return name[:34]


_ORG_RANKS = {"subsp", "subsp.", "subspecies", "ssp", "ssp.", "var", "var.",
              "sp", "sp.", "str", "str.", "substr", "substr.", "strain",
              "biovar", "serovar", "pv", "pv.", "f.", "form"}
_ORG_STRAINY = ("str", "substr", "strain", "biovar", "serovar", "pv")


def italic_organism(org: str) -> str:
    """Format an organism name for a matplotlib title with the scientific name in
    ITALICS (genus, species, subspecies epithet, and 'Candidatus') and strain
    designations upright, with the assembly accession in parentheses if present.
    Returns a mathtext string: the italic parts are '$\\mathit{...}$' (rendered
    non-bold via the custom mathtext fontset); strain/accession are plain
    (upright) text so hyphens/digits render normally. Non-standard names (e.g.
    'xyz genome') degrade gracefully -- the leading word(s) may be italicised,
    which is harmless for a reference figure (captions carry the real name)."""
    acc = re.search(r"(GC[AF]_\d+\.\d+)", org or "")
    name = re.sub(r"_GC[AF]_.*$", "", org or "").replace("_", " ").strip()
    toks = name.split()
    if not toks:
        return (org or "") + (f" ({acc.group(1)})" if acc else "")
    styles, strain = [], False
    for k, t in enumerate(toks):
        tl = t.lower()
        if tl in _ORG_RANKS:
            styles.append("rm")
            if tl.startswith(_ORG_STRAINY):
                strain = True
            continue
        if strain:
            styles.append("rm")
            continue
        if k == 0 or tl == "candidatus":
            styles.append("it")
            continue
        if re.fullmatch(r"[a-z-]+", t):   # lowercase epithet (original case)
            styles.append("it")
            continue
        styles.append("rm")               # uppercase/digit token = strain designation
        strain = True
    parts, i = [], 0
    while i < len(toks):
        if styles[i] == "it":
            j = i
            while j < len(toks) and styles[j] == "it":
                j += 1
            parts.append(r"$\mathit{" + r"\ ".join(toks[i:j]) + "}$")
            i = j
        else:
            parts.append(toks[i])
            i += 1
    label = " ".join(parts)
    if acc:
        label += f" ({acc.group(1)})"
    return label


def organism_segments(org: str):
    """Split an organism name into [(text, is_italic)] segments: scientific name
    parts italic, strain/accession upright. Used by draw_organism_line() to
    render real (non-bold) italics via offsetbox instead of mathtext."""
    acc = re.search(r"(GC[AF]_\d+\.\d+)", org or "")
    name = re.sub(r"_GC[AF]_.*$", "", org or "").replace("_", " ").strip()
    toks = name.split()
    segs = []
    if toks:
        styles, strain = [], False
        for k, t in enumerate(toks):
            tl = t.lower()
            if tl in _ORG_RANKS:
                styles.append(False)
                if tl.startswith(_ORG_STRAINY):
                    strain = True
            elif strain:
                styles.append(False)
            elif k == 0 or tl == "candidatus":
                styles.append(True)
            elif re.fullmatch(r"[a-z-]+", t):
                styles.append(True)
            else:
                styles.append(False); strain = True
        i = 0
        while i < len(toks):
            st = styles[i]; j = i
            while j < len(toks) and styles[j] == st:
                j += 1
            segs.append((" ".join(toks[i:j]), st))
            i = j
    if acc:
        segs.append((f"({acc.group(1)})", False))
    if not segs:
        segs = [(org or "", False)]
    return segs


def _organism_hpacker(org: str, fontsize: float):
    from matplotlib.offsetbox import TextArea, HPacker
    segs = organism_segments(org)
    boxes = [TextArea(t, textprops=dict(
                fontstyle="italic" if it else "normal", fontweight="normal",
                fontsize=fontsize, color=INK)) for t, it in segs if t]
    return HPacker(children=boxes, align="baseline", pad=0, sep=fontsize * 0.30)


def draw_organism_line(fig, org: str, y: float = 0.945, fontsize: float = 15) -> None:
    """Draw the organism name centred at figure-fraction y as a plain Text
    artist (robust under tight-bbox cropping, unlike an offsetbox): the
    scientific name is italicised via mathtext and the whole line is non-bold."""
    t = fig.text(0.5, y, italic_organism(org), ha="center", va="top",
                 fontsize=fontsize, fontweight="normal", color=INK)
    fig._report_extra_artists = getattr(fig, "_report_extra_artists", []) + [t]


def draw_title(fig, description: str, org: str, fontsize_desc: float = 16,
               fontsize_org: float = 13.5, y: float = 0.998) -> None:
    """Stack a bold description over the (non-bold, italic scientific name)
    organism line as ONE top-centred offsetbox -- so the two never overlap or
    clip regardless of figure height / tight-bbox cropping."""
    from matplotlib.offsetbox import TextArea, VPacker, AnnotationBbox
    desc = TextArea(description, textprops=dict(fontweight="bold",
                    fontsize=fontsize_desc, color=INK))
    stack = VPacker(children=[desc, _organism_hpacker(org, fontsize_org)],
                    align="center", pad=0, sep=fontsize_org * 0.5)
    ab = AnnotationBbox(stack, (0.5, y), xycoords="figure fraction", frameon=False,
                        box_alignment=(0.5, 1.0))
    fig.add_artist(ab)
    fig._report_extra_artists = getattr(fig, "_report_extra_artists", []) + [ab]


# ---------------------------------------------------------------------------
# Data loading -- FINISHED scoring outputs only.
# ---------------------------------------------------------------------------
_CONF_FINAL = "scoring/scored-labeled-genes-confidence-final.tsv"
_LABELED = "labeling/labeled-genes.tsv"
_NOT_IN_OPERON = "NOT_IN_AN_OPERON"


def is_operon(oid) -> bool:
    """True only for a real operon id. Real operons all start with 'operon_';
    NOT_IN_AN_OPERON, NOT_APPLICABLE_NON_CODING, blanks, etc. are NOT operons."""
    return bool(oid) and str(oid).startswith("operon_")

_NUMERIC_COLS = [
    "c1_score", "c2_score_from_operon_probability", "c3_score", "c4_score",
    "preliminary_confidence_c1_c4", "final_confidence_operon_context",
    "confidence_score", "operon_member_count", "operon_gene_position_in_operon",
]
_GENE_ID_RE = re.compile(r"^(.*)_(\d+)([+-])(\d+)$")


def parse_contig(gene_id: str) -> str:
    m = _GENE_ID_RE.match(gene_id or "")
    return m.group(1) if m else (gene_id or "")


def discover_organisms(run_root: Path) -> list[str]:
    """Organism stems under the run that have a confidence-final scoring file."""
    run_root = Path(run_root)
    out = []
    for p in sorted(run_root.glob("*/" + _CONF_FINAL)):
        out.append(p.parents[1].name)
    return out


def _coerce_numeric(df: pd.DataFrame) -> pd.DataFrame:
    for c in _NUMERIC_COLS:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def load_organism_genes(run_root: Path, organism: str) -> pd.DataFrame:
    """One row per gene for `organism`: confidence-final scores + operon info,
    joined to genomic coordinates (start/end/strand/contig) from labeled-genes.
    Adds `clean_desc`, `uninformative`, `in_operon`, `contig`."""
    run_root = Path(run_root)
    base = run_root / organism
    conf = pd.read_csv(base / _CONF_FINAL, sep="\t", dtype=str,
                       keep_default_na=False, engine="python")
    conf = _coerce_numeric(conf)

    coords_cols = ["feature_id", "gene_id", "gene_start", "gene_end", "RAST_strand"]
    lab = pd.read_csv(base / _LABELED, sep="\t", dtype=str,
                      keep_default_na=False, engine="python",
                      usecols=lambda c: c in coords_cols)
    for c in ("gene_start", "gene_end"):
        if c in lab.columns:
            lab[c] = pd.to_numeric(lab[c], errors="coerce")
    df = conf.merge(lab, on="feature_id", how="left")

    df["contig"] = df.get("gene_id", "").map(parse_contig)
    df["clean_desc"] = df["best_consensus_product_descriptor"].map(clean_descriptor)
    df["uninformative"] = df["best_consensus_product_descriptor"].map(is_uninformative)
    # A gene is IN an operon only if its operon_id is a real operon id (they all
    # start with "operon_"). Everything else -- NOT_IN_AN_OPERON,
    # NOT_APPLICABLE_NON_CODING (non-coding genes), blanks -- is NOT an operon.
    oid = df.get("operon_id", pd.Series([""] * len(df))).astype(str)
    df["in_operon"] = oid.str.startswith("operon_")
    df["organism"] = organism
    return df


def load_all_genes(run_root: Path, organisms: list[str] | None = None) -> pd.DataFrame:
    run_root = Path(run_root)
    organisms = organisms or discover_organisms(run_root)
    frames = []
    for org in organisms:
        try:
            frames.append(load_organism_genes(run_root, org))
        except Exception as e:  # pragma: no cover
            print(f"[reportfig] skip {org}: {e}", file=sys.stderr)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


# ---- operons with ordered members + coordinates -------------------------------
def build_operons(genes: pd.DataFrame) -> pd.DataFrame:
    """One row per operon (in the given genes frame): ordered members, member
    labels, coordinates, size, and the `members_in_order` join key used by the
    depot operon-fingerprint database ("label1 -> label2 -> ...")."""
    ops = genes[genes["in_operon"]].copy()
    if ops.empty:
        return pd.DataFrame(columns=["organism", "operon_id", "members_in_order",
                                     "size", "member_labels", "feature_ids"])
    ops["_pos"] = pd.to_numeric(ops.get("operon_gene_position_in_operon"),
                                errors="coerce").fillna(0)
    rows = []
    for (org, oid), g in ops.groupby(["organism", "operon_id"]):
        g = g.sort_values("_pos")
        labels = g["best_consensus_product_descriptor"].tolist()
        rows.append({
            "organism": org,
            "operon_id": oid,
            "members_in_order": " -> ".join(labels),
            "size": len(labels),
            "member_labels": labels,
            "clean_labels": g["clean_desc"].tolist(),
            "sources": g.get("product_descriptor_source",
                             pd.Series([""] * len(g))).tolist(),
            "needs_reviews": g.get("needs_review",
                                   pd.Series([""] * len(g))).tolist(),
            "feature_types": g.get("feature_type",
                                   pd.Series([""] * len(g))).tolist(),
            "feature_ids": g["feature_id"].tolist(),
            "starts": g["gene_start"].tolist(),
            "ends": g["gene_end"].tolist(),
            "strands": g.get("RAST_strand", pd.Series([""] * len(g))).tolist(),
            "contigs": g["contig"].tolist(),
            "confidences": pd.to_numeric(
                g.get("confidence_score", pd.Series([np.nan] * len(g))),
                errors="coerce").tolist(),
            "preliminaries": pd.to_numeric(
                g.get("preliminary_confidence_c1_c4", pd.Series([np.nan] * len(g))),
                errors="coerce").tolist(),
            "operon_adjusteds": pd.to_numeric(
                g.get("final_confidence_operon_context", pd.Series([np.nan] * len(g))),
                errors="coerce").tolist(),
            "operon_adjusteds_hybrid": pd.to_numeric(
                g.get("final_confidence_operon_context_hybrid", pd.Series([np.nan] * len(g))),
                errors="coerce").tolist(),
            "c3s": pd.to_numeric(
                g.get("c3_score", pd.Series([np.nan] * len(g))),
                errors="coerce").tolist(),
            "c3s_hybrid": pd.to_numeric(
                g.get("c3_score_operon_context_hybrid", pd.Series([np.nan] * len(g))),
                errors="coerce").tolist(),
            "c2s": pd.to_numeric(
                g.get("c2_score_from_operon_probability", pd.Series([np.nan] * len(g))),
                errors="coerce").tolist(),
            "c1s": pd.to_numeric(
                g.get("c1_score", pd.Series([np.nan] * len(g))),
                errors="coerce").tolist(),
            "c4s": pd.to_numeric(
                g.get("c4_score", pd.Series([np.nan] * len(g))),
                errors="coerce").tolist(),
            "review_reasons": g.get("needs_review_reason",
                                    pd.Series([""] * len(g))).tolist(),
        })
    return pd.DataFrame(rows)


# ---- depot pangenome operon recurrence ----------------------------------------
def load_operon_recurrence(depot_db: Path, restrict_to=None) -> dict[str, dict]:
    """members_in_order -> {label_frequency, organism_count, organisms} from the
    depot operon-fingerprint label-ordered database (the pangenome operon pool).
    Read-only.

    `restrict_to`: an optional iterable of organism names. When given, recurrence
    is SCOPED to only those organisms -- each operon's `organisms`/`organism_count`
    is intersected with the set, so "in K pangenome genomes" counts only THIS
    project's organisms (e.g. the genomes in input-user/) and never the extra
    reference genomes that also accumulate in the shared depot DB. Dynamic: the
    caller reads the set fresh each run, so it grows automatically as more genomes
    are added. An empty/None set means no scoping (the whole DB)."""
    depot_db = Path(depot_db)
    keep = set(restrict_to) if restrict_to else None
    out: dict[str, dict] = {}
    if not depot_db.is_file():
        print(f"[reportfig] WARNING: operon recurrence DB not found: {depot_db}",
              file=sys.stderr)
        return out
    df = pd.read_csv(depot_db, sep="\t", dtype=str, keep_default_na=False,
                     engine="python")
    for _, r in df.iterrows():
        mio = r.get("members_in_order", "")
        # DEDUP the organism list: recurrence counts DISTINCT organisms, so an
        # organism can never be double-counted even if the depot OCC DB were ever
        # to list the same organism twice for an operon (guards against a
        # re-processing bug inflating "in K pangenome genomes").
        orgs = sorted({x for x in (r.get("organisms", "") or "").split("|") if x})
        if keep is not None:
            orgs = [o for o in orgs if o in keep]
        n = len(orgs)
        rec = {
            # when scoped, the DB-wide occurrence count is not per-organism
            # decomposable, so occurrences collapse to the scoped organism count
            "label_frequency": (n if keep is not None
                                else int(r.get("fingerprint_label_frequency", 0) or 0)),
            "organism_count": n,
            "organisms": orgs,
        }
        # keep the richest record if the same text appears under >1 hash
        prev = out.get(mio)
        if prev is None or rec["organism_count"] > prev["organism_count"]:
            out[mio] = rec
    return out


# Pangenome recurrence is scoped to the organisms discovered in the RUN'S output
# folder (discover_organisms) -- the organisms actually scored in that run and the
# exact OCC pool its C3 was computed against. This is independent of where the
# user keeps their input genomes, so it works for any input path.


DEFAULT_OPERON_DB = Path(
    "/depot/lindems/data/margie/fingerprint-database/"
    "operon-fingerprint-database-label-ordered.tsv"
)
DEFAULT_OCC_REFERENCE = Path(
    "/depot/lindems/data/margie/operon-database/occ_reference.pkl"
)


def load_occ_organisms(occ_reference=DEFAULT_OCC_REFERENCE):
    """The organisms in the OCC reference (occ_reference.pkl's `organisms_added`)
    -- i.e. the ACTUAL pool the C3 scores were computed against, the persistent
    baseline that accumulates across runs. This is the correct thing to scope the
    figures' recurrence + pool caption to: unlike discover_organisms(run) (which
    only sees the genomes scored so far in ONE run, so mid-run it under-reports,
    e.g. "4 genomes" while C3 actually used all 21), this reflects the real OCC.
    Returns a set of organism names, or None if unreadable (callers fall back to
    the run's organisms)."""
    try:
        import pickle
        with open(occ_reference, "rb") as fh:
            ref = pickle.load(fh)
        orgs = {str(o) for o in (ref.get("organisms_added") or []) if o}
        return orgs or None
    except Exception:
        return None


_POOL_STAT_COLS = ("total_genes", "operonic_genes", "singleton_genes",
                   "n_operons", "n_informative_operons", "n_uninformative_operons")


def load_pool_stats(occ_reference=DEFAULT_OCC_REFERENCE):
    """Per-genome pool stats from the OCC's .genome_stats.tsv sidecar:
    {organism: {total_genes, operonic_genes, singleton_genes, n_operons,
    n_informative_operons, n_uninformative_operons}}. Read straight from the
    sidecar TSV (no run-folder dependency), so it is complete even mid-run.
    Empty dict if the sidecar is absent/unreadable."""
    import csv
    path = Path(str(occ_reference) + ".genome_stats.tsv")
    out = {}
    try:
        with open(path, newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                org = (row.get("organism") or "").strip()
                if org:
                    out[org] = {c: int(float(row.get(c, 0) or 0))
                                for c in _POOL_STAT_COLS}
    except Exception:
        return {}
    return out


def aggregate_pool_stats(stats_map, organisms):
    """Sum per-genome pool stats over ``organisms`` (those present in stats_map).
    Returns a dict over the stat columns plus ``n_genomes`` (how many had stats)."""
    agg = {c: 0 for c in _POOL_STAT_COLS}
    n = 0
    for org in organisms:
        s = stats_map.get(org)
        if s:
            n += 1
            for c in _POOL_STAT_COLS:
                agg[c] += s[c]
    agg["n_genomes"] = n
    return agg


def rel_or_host(path, run_root) -> str:
    """Display a source path for the on-figure footer: relative to the run's
    output directory (tagged [output]) when it lives inside the run, else the
    absolute host path (tagged [host]) for compute-host reference files."""
    p, run_root = Path(path), Path(run_root)
    try:
        return f"{p.relative_to(run_root)}  [output]"
    except ValueError:
        return f"{p}  [host]"


def organism_source_lines(run_root, organism: str, coords: bool = False,
                          operon_db=None) -> list[str]:
    """Footer source list for a per-organism figure."""
    run_root = Path(run_root)
    lines = [rel_or_host(run_root / organism / _CONF_FINAL, run_root)]
    if coords:
        lines.append(rel_or_host(run_root / organism / _LABELED, run_root))
    if operon_db is not None:
        lines.append(rel_or_host(operon_db, run_root))
    return lines


def global_source_lines(run_root, n_org: int, scored: bool = True,
                        operon_db=None) -> list[str]:
    """Footer source list for a run-level (pangenome) figure."""
    lines = []
    if scored:
        fname = _CONF_FINAL.split("/")[-1]
        lines.append(f"*/scoring/{fname}  (all {n_org} genomes)  [output]")
    if operon_db is not None:
        lines.append(rel_or_host(operon_db, run_root))
    return lines


# Compact per-column glossary, rendered LEFT-aligned (draw_method_note) so the
# definitions read as a clean two/three-line key rather than a centred paragraph.
# The "\n" splits are honoured verbatim (draw_method_note wraps only over-long
# lines), so each column's term = definition stays grouped.
OPERON_CORRECTION_NOTE = (
    "Columns —  "
    "C1 = tool coverage/agreement  |  "
    "C2 = pairwise operon probability with the adjacent gene (UniOP)  |  "
    "C3(adj) = geomean(ρ) over the operon's adjacent informative pairs (operon-level; all members share)  |  "
    "C3(hyb) = max(best adjacent ρ, best co-member ρ), per gene  |  "
    "C4 = EC agreement (< 1 = conflict)\n"
    "prelim = C1×C4  |  "
    "operon boost = C2×C3, boost-only, shown adj/hyb (raises an uncertain gene, never lowers; novel C3→0 has been kept neutral and can only be penalized if the operon database is sufficient: hence, left for future work)  |  "
    "final = clip(prelim + boost, 0, 1), shown adj/hyb  |  "
    "review? = manual-check flag  |  "
    "review reason = EC conflict / ambiguous operon / low confidence (final < 0.5)  |  "
    "ρ = cross-genome pair reliability (Jeffreys posterior × enrichment)"
)

OPERON_PENALTY_NOTE = (
    "Reading the trail: ‘operon boost’ here is NEGATIVE — a penalty = C2·conflict, where conflict is a "
    "DESCRIPTOR contradiction: across the pooled genomes this gene’s operon slot carries a DIFFERENT "
    "functional descriptor by consensus, so the functional call we assigned is the contradicted minority. "
    "C2 = operon probability only gates it (a doubtful operon shrinks the penalty; it can never create one). "
    "Crucially this is NOT a penalty for novelty: a merely-unseen operon (C3 → 0, no conflict) is neutral. "
    "Only positive cross-genome contradiction of the descriptor lowers the score — evidence the functional "
    "annotation, not the gene’s placement, is likely wrong. ‘review?’ reads ‘yes’ only when the drop is "
    "material (≥ 0.1)."
)


_STRW_CACHE: dict = {}   # (word, fontsize) -> width in figure-fraction units


def _fig_strwidth(fig, renderer, s, fontsize):
    """Width of string `s` (at `fontsize`) as a fraction of the figure width.
    Cached: the glossary text is identical on every atlas page, so after the
    first page every measurement is a dict hit (keeps the atlas fast)."""
    key = (s, fontsize)
    w = _STRW_CACHE.get(key)
    if w is not None:
        return w
    t = fig.text(0, 0, s, fontsize=fontsize)
    try:
        bb = t.get_window_extent(renderer=renderer)
        w = bb.width / (fig.get_figwidth() * fig.dpi)
    except Exception:
        w = len(s) * 0.0128 * fontsize / fig.get_figwidth()
    t.remove()
    _STRW_CACHE[key] = w
    return w


def draw_method_note(fig, text: str, fontsize: float = 8.2,
                     left: float = 0.015, right: float = 0.985) -> None:
    """Print a methodology / column-glossary note just below the plot content
    (above the sources footer). Rendered FULLY JUSTIFIED — each line's words are
    placed individually so the block has clean left AND right edges (last line of
    each ``\\n`` paragraph is left-aligned). Measured placement + registered so
    the tight crop keeps it."""
    if not text:
        return
    fh = max(fig.get_figheight(), 3.0)
    try:
        fig.draw_without_rendering()
        renderer = fig.canvas.get_renderer()
        extras = getattr(fig, "_report_extra_artists", None)
        low_in = fig.get_tightbbox(renderer, bbox_extra_artists=extras).y0
    except Exception:
        renderer, low_in = None, 0.25

    def sw(s):
        return _fig_strwidth(fig, renderer, s, fontsize)

    avail = right - left
    space_w = sw(" ")
    line_h = (fontsize * 1.55 / 72.0) / fh          # line pitch in figure fraction
    y = (low_in - 0.22) / fh
    artists = []
    for para in text.split("\n"):
        words = para.split()
        if not words:
            y -= line_h
            continue
        # greedy word-wrap to the available width
        rows, cur, cur_w = [], [], 0.0
        for wd in words:
            ww = sw(wd)
            add = ww + (space_w if cur else 0.0)
            if cur and cur_w + add > avail:
                rows.append(cur)
                cur, cur_w = [wd], ww
            else:
                cur.append(wd)
                cur_w += add
        if cur:
            rows.append(cur)
        for ri, row in enumerate(rows):
            last = (ri == len(rows) - 1)
            widths = [sw(wd) for wd in row]
            if last or len(row) == 1:
                gap = space_w                        # ragged last line: normal spacing
            else:
                gap = (avail - sum(widths)) / (len(row) - 1)   # justify: spread slack
            x = left
            for wd, ww in zip(row, widths):
                artists.append(fig.text(x, y, wd, ha="left", va="top",
                                        fontsize=fontsize, color=_DETAIL_PURPLE))
                x += ww + gap
            y -= line_h                              # advance one line per wrapped row
    fig._report_extra_artists = getattr(fig, "_report_extra_artists", []) + artists


def draw_sources_footer(fig, run_root, sources: list[str], fontsize: float = 6.8) -> None:
    """Print a small, non-distracting footer at the bottom-centre of the figure
    naming the exact files that fed this graph, so a viewer can audit them.
    Paths are shown relative to the run's output directory ([output]) or as the
    absolute compute-host path ([host]). Reserves a bottom band so it never
    overlaps the plot, and registers itself so the tight crop keeps it."""
    if not sources:
        return
    run_root = Path(run_root)
    # The label ends with ":" -- the ONLY ":" used as a separator in these figures;
    # everything else (here and elsewhere) separates with "|".
    header = (f"Files used to build this figure:  paths relative to output/{run_root.name}/  "
              f"|  [host] = reference file on the compute host")
    fw, fh = fig.get_figwidth(), max(fig.get_figheight(), 3.0)
    approx = max(60, int(fw / (0.011 * fontsize)))
    body = "   |   ".join(sources)
    lines = [header] + (textwrap.wrap(body, width=approx) or [body])
    # Measure the current lowest content (x tick labels / axis title / table) and
    # drop the footer just below it. No subplots_adjust -- the explicit-bbox crop
    # in savefig() extends downward to include the footer, so it can even sit
    # below y=0. This is robust to multi-line tick labels and any figure height.
    try:
        fig.draw_without_rendering()
        r = fig.canvas.get_renderer()
        extras = getattr(fig, "_report_extra_artists", None)
        low_in = fig.get_tightbbox(r, bbox_extra_artists=extras).y0
    except Exception:
        low_in = 0.25
    y_top = (low_in - 0.18) / fh                 # 0.18" gap below the lowest content
    t = fig.text(0.5, y_top, "\n".join(lines), ha="center", va="top",
                 fontsize=fontsize, color=_MAGENTA, fontweight="normal",
                 linespacing=1.4)
    fig._report_extra_artists = getattr(fig, "_report_extra_artists", []) + [t]


def _describe_figure_tsv(name: str) -> str:
    fig = name.split("_", 1)[0]
    rest = name.split("_", 1)[1].rsplit(".", 1)[0].replace("_", " ") if "_" in name else ""
    return f"exact data table plotted in {fig}: {rest}".rstrip(": ")


def write_sources_manifest(outdir: Path, run_root: Path, organisms: list[str],
                           operon_db: Path,
                           occ_reference: Path = DEFAULT_OCC_REFERENCE) -> None:
    """Write `figure-sources.tsv` into the figures folder: every ORIGINAL file
    consulted to build these figures, with its absolute location, whether it was
    found, its size/timestamp, and what it provided -- plus every processed
    fig*_*.tsv companion (the exact numbers plotted). Lets a reviewer trace each
    figure back to the scored file, the depot databases and the OCC reference."""
    from datetime import datetime
    outdir, run_root = Path(outdir), Path(run_root)
    rows = []

    def add(role: str, path, description: str) -> None:
        p = Path(path)
        try:
            ok = p.exists()
            size = p.stat().st_size if ok else ""
            mtime = (datetime.fromtimestamp(p.stat().st_mtime).isoformat(timespec="seconds")
                     if ok else "")
        except OSError:
            ok, size, mtime = False, "", ""
        rows.append({"role": role, "item": p.name, "path": str(p),
                     "found": "yes" if ok else "MISSING",
                     "size_bytes": size, "modified": mtime,
                     "description": description})

    for org in organisms:
        base = run_root / org
        add("input: scored genes", base / _CONF_FINAL,
            f"[{org}] final per-gene confidence table (C1-C4, preliminary, "
            "operon-adjusted, final score, tier) with operon_id and member "
            "counts -- the values shown in the plots")
        add("input: gene coordinates", base / _LABELED,
            f"[{org}] gene genomic coordinates (start/end/strand) and product "
            "descriptors, joined in for the gene-arrow maps and locations")
    add("database: operon recurrence", operon_db,
        "depot pangenome operon-fingerprint database (label-ordered): how many "
        "genomes each operon recurs across and how often -- drives the "
        "reproduced / unique / most-conserved operon figures")
    add("reference: operon co-occurrence (OCC)", occ_reference,
        "depot operon co-occurrence reference behind the operon-probability (C2) "
        "and operon-context (C3) scores that scoring wrote into the scored file")
    for tsv in sorted(outdir.glob("fig*_*.tsv")):
        add("figure data (plotted)", tsv, _describe_figure_tsv(tsv.name))

    df = pd.DataFrame(rows, columns=["role", "item", "path", "found",
                                     "size_bytes", "modified", "description"])
    write_tsv(df, outdir / "figure-sources.tsv")


# ---- small numeric helpers ----------------------------------------------------
def pca_svd(X: np.ndarray):
    """Standardized PCA via SVD (no sklearn). Returns (scores, loadings,
    explained_variance_ratio). Columns of X are variables."""
    Xc = X - np.nanmean(X, axis=0)
    sd = np.nanstd(Xc, axis=0)
    sd[sd == 0] = 1.0
    Xs = Xc / sd
    Xs = np.nan_to_num(Xs)
    U, S, Vt = np.linalg.svd(Xs, full_matrices=False)
    scores = U * S
    evr = (S ** 2) / np.sum(S ** 2)
    return scores, Vt.T, evr


def size_bin(n) -> str:
    """Human-readable operon-size bands (member counts)."""
    try:
        n = int(n)
    except (TypeError, ValueError):
        return "unknown"
    if n <= 1:
        return "1 (not in operon)"
    if n == 2:
        return "2"
    if n <= 4:
        return "3-4"
    if n <= 8:
        return "5-8"
    if n <= 20:
        return "9-20"
    return "21+"


SIZE_BIN_ORDER = ["2", "3-4", "5-8", "9-20", "21+"]


def box_by_bin(ax, df, cat_col, value_col, order, color, ylabel, xlabel):
    """Boxplot of value_col grouped by an ordered categorical column, with the
    MEAN of each bin overlaid as a connected line so mean vs median can be read
    off the same axis (both share the identical 0-1 scale). The box shows the
    median (centre line) + interquartile spread; a skew between the two is the
    point. Returns (labels, ns, medians, means) for the companion TSV."""
    from matplotlib.lines import Line2D
    data, labels, ns = [], [], []
    for b in order:
        vals = pd.to_numeric(df.loc[df[cat_col] == b, value_col],
                             errors="coerce").dropna()
        if len(vals):
            data.append(vals.values); labels.append(b); ns.append(len(vals))
    if not data:
        ax.axis("off"); return [], [], [], []
    xs = list(range(len(data)))
    bp = ax.boxplot(data, positions=xs, widths=0.6,
                    patch_artist=True, showfliers=False,
                    medianprops=dict(color=INK, lw=2.0),
                    whiskerprops=dict(color="#555555"), capprops=dict(color="#555555"))
    # each bin gets its own vibrant colour
    for i, patch in enumerate(bp["boxes"]):
        patch.set_facecolor(CATEGORICAL[i % len(CATEGORICAL)])
        patch.set_alpha(0.85); patch.set_edgecolor("white"); patch.set_linewidth(1.5)
    medians = [float(np.median(x)) for x in data]
    means = [float(np.mean(x)) for x in data]
    # mean overlay: black line + white diamonds (reads clearly on top of boxes)
    ax.plot(xs, means, color=INK, lw=1.6, marker="D", markersize=8, zorder=6,
            markerfacecolor="white", markeredgecolor=INK, markeredgewidth=1.6)
    ax.set_xticks(xs)
    ax.set_xticklabels([f"{l}\n(n={n:,})" for l, n in zip(labels, ns)], fontsize=9)
    ax.set_ylabel(ylabel); ax.set_xlabel(xlabel)
    ax.grid(False)
    handles = [Line2D([0], [0], color=INK, lw=2.0, label="median (box centre)"),
               Line2D([0], [0], color=INK, lw=1.6, marker="D", markersize=8,
                      markerfacecolor="white", markeredgecolor=INK, label="mean")]
    ax.legend(handles=handles, loc="lower right", fontsize=8.5, frameon=True,
              framealpha=0.92, edgecolor="#cccccc")
    return labels, ns, medians, means


# ---------------------------------------------------------------------------
# Gene-track drawing (the "real operon neighbourhood" idiom) + descriptor index.
# Genes are drawn as coloured nodes tagged A, B, C, …; the FULL descriptor for
# each tag (with its final confidence score) is listed in an index -- so gene
# names are never truncated with "…".
# ---------------------------------------------------------------------------
def tag_letter(i: int) -> str:
    """Gene tag for the operon diagrams: 0->'1', 1->'2', … Numbered (not lettered)
    so it stays clear and scales past 26 -- operons run up to ~79 genes, where
    A/B/…/Z/AA/AB gets unreadable. The 'gene' table column shows the same number."""
    return str(i + 1)


# Target PHYSICAL size of a gene glyph, in inches. draw_gene_track converts
# these to data units from the axis's real width/height so every arrow -- in a
# dense fig05 track or a sparse 2-gene gallery row -- looks the SAME: a compact
# block arrow with a short (non-"rocket") head, never a thin spear or a fat wedge.
_ARROW_W_IN = 0.74     # total glyph LENGTH (horizontal) -- shorter/stubbier, not stretched
_ARROW_HLEN_IN = 0.16  # head length -- SHORT pointer (stubby head, long body)
_ARROW_BODY_IN = 0.18  # body thickness (vertical) -- fixed PHYSICAL inches
_ARROW_HEAD_IN = 0.36  # head thickness (vertical, wider than body) -- fixed inches
_TRACK_YRANGE = 1.9    # data height ONE arrow row maps to (== row_pitch below).
_TRACK_ROW_IN = 1.28   # physical inches allotted to ONE arrow row (atlas contract;
                       # the atlas driver uses this same value for its row height).
# Inches per data-y-unit that EVERY arrow track is pinned to. The arrow glyph is
# sized in fixed data units from this, and pin_track_scale() resets each track's
# ylim AFTER figure layout so its real inches-per-data-unit equals this exactly --
# so a 1-row and a 6-row operon get identically thick arrows on every page,
# immune to tight_layout's per-page rescaling.
_TRACK_IN_PER_Y = _TRACK_ROW_IN / _TRACK_YRANGE


def _arrow_geometry(ax, span: int, nrows: int = 1):
    """Data-unit glyph geometry (hw, body_ht, head_ht, head_len). The VERTICAL
    thickness is a fixed data size derived from the pinned _TRACK_IN_PER_Y (NOT
    the measured axis height, which tight_layout distorts differently per page);
    pin_track_scale() then makes the rendered inches-per-data-unit match, so the
    arrow is exactly _ARROW_HEAD_IN/_ARROW_BODY_IN inches tall on every page. The
    HORIZONTAL size still uses the measured axis width (full 18in, stable)."""
    fig = ax.figure
    try:
        pos = ax.get_position()
        ax_w_in = max(pos.width * fig.get_figwidth(), 1.0)
    except Exception:
        ax_w_in = 12.0
    slot_in = ax_w_in / max(span, 1)
    in_per_y = _TRACK_IN_PER_Y
    hw = min(0.47, (_ARROW_W_IN / slot_in) / 2.0)
    hlen = min(hw * 0.8, _ARROW_HLEN_IN / slot_in)
    ht = min(0.42, (_ARROW_HEAD_IN / 2.0) / in_per_y)
    bt = min(ht * 0.6, (_ARROW_BODY_IN / 2.0) / in_per_y)
    return hw, bt, ht, hlen


def pin_track_scale(ax) -> None:
    """Re-pin an arrow track's ylim AFTER figure layout so its rendered
    inches-per-data-unit equals _TRACK_IN_PER_Y exactly. draw_gene_track sizes
    arrows in fixed data units for that scale, but tight_layout/subplots_adjust
    resize each page's axes by a different factor -- so without this, wrapped
    (multi-row) pages still render arrows slightly thinner than single-row pages.
    Call once per track axis after the figure's final layout is set."""
    info = getattr(ax, "_track_scale", None)
    if not info:
        return
    fig = ax.figure
    try:
        ax_h_in = ax.get_position().height * fig.get_figheight()
    except Exception:
        return
    if ax_h_in <= 0:
        return
    rng = ax_h_in / _TRACK_IN_PER_Y
    # never frame tighter than the drawn content (guards against clipping a
    # wrapped row if layout squeezed the axis below its nominal height)
    rng = max(rng, info["nrows"] * _TRACK_YRANGE)
    ax.set_ylim(info["y_top"] - rng, info["y_top"])


def member_descriptor(m) -> str:
    """Final table descriptor for a gene-track member: the cleaned product
    descriptor prefixed with its annotation source (PGAP/UNIPROT/…) when known.
    Used both to DRAW the table and to SIZE it, so they never disagree."""
    full = clean_descriptor(m.get("label", "")) or (m.get("label", "") or "—")
    src = (m.get("source") or "").strip()
    if src and full != "—":
        full = f"{src}: {full}"
    return full


def member_table_units(members, desc_wrap: int = None) -> float:
    """Wrapped-line count a gene table needs for these members (2-line header +
    each member's wrapped, SOURCE-PREFIXED descriptor) -- size the table axis
    with this so long 'PGAP: …' names never overrun their row."""
    w = desc_wrap or _TABLE_DESC_WRAP
    return _TABLE_HEADER_UNITS + sum(len(wrap_desc(member_descriptor(m), w)) for m in members)


def _gene_arrow(cx, y, strand, col, geom):
    """Classic block-arrow gene glyph pointing 5'→3' (right for +, left for −):
    a rectangular body with a wider triangular head that tapers to a point."""
    hw, bt, ht, hl = geom
    if str(strand) == "-":
        # tip at left
        tip, sh = cx - hw, cx - hw + hl
        pts = [(cx + hw, y + bt), (sh, y + bt), (sh, y + ht), (tip, y),
               (sh, y - ht), (sh, y - bt), (cx + hw, y - bt)]
    else:
        # tip at right
        tip, sh = cx + hw, cx + hw - hl
        pts = [(cx - hw, y + bt), (sh, y + bt), (sh, y + ht), (tip, y),
               (sh, y - ht), (sh, y - bt), (cx - hw, y - bt)]
    return Polygon(pts, closed=True, facecolor=col, edgecolor="#222222",
                   linewidth=1.1, joinstyle="miter", zorder=3)


def draw_gene_track(ax, members: list[dict], badge_per_operon: dict | None = None,
                    show_gaps: bool = True, letter_offset: int = 0,
                    node_size: int = 760, strand_rows: bool = False,
                    min_span: int | None = None, per_row: int | None = None,
                    tag_fs: float = 12.0, gap_fs: float = 11.0):
    """Draw genes as a classical single-line gene-arrow map (gggenes/clinker
    style): one horizontal genome backbone, each gene a block arrow pointing
    5'→3' along its own coding strand, coloured per gene and tagged A/B/C; a
    light rounded band groups the genes of one operon; the intergenic gap in bp
    sits under the backbone between neighbours. `min_span` fixes the number of
    gene slots the axis spans (genes are centred, backbone runs full width as
    flanking genome) so arrows keep a consistent, compact size regardless of how
    many genes a row has -- otherwise a 2-gene row stretched across a wide axis
    would draw absurdly long arrows. `strand_rows` is accepted for
    backward-compatibility but ignored (strand is shown by arrow direction, so a
    single line is used). `per_row` wraps a long operon onto multiple stacked
    rows of at most that many genes each (arrows keep a consistent size and tags
    A/B/C… run on across rows) so operons far larger than one line -- up to 79
    genes -- render losslessly; None keeps the classic single row. Returns
    5-tuples (tag, descriptor, score, location, colour) for render_gene_table()."""
    n = len(members)
    if n == 0:
        ax.axis("off")
        return []
    gene_col = [CATEGORICAL[(letter_offset + i) % len(CATEGORICAL)] for i in range(n)]
    # row layout: wrap onto rows of at most `per_row` genes when requested and
    # the operon is longer than one row; otherwise everything sits on one row
    # (per == n, nrows == 1) which reproduces the classic single-line map exactly.
    per = per_row if (per_row and n > per_row) else n
    per = max(int(per), 1)
    nrows = (n + per - 1) // per
    span = max(per, min_span) if min_span else per
    off = (span - per) / 2.0                  # centre the per-column grid
    x0, x1 = -0.60, span - 0.40

    geom = _arrow_geometry(ax, span, nrows)
    ht = geom[2]                              # head half-thickness (data units)
    band_h = ht + 0.16                        # operon band half-height
    tag_y = ht + 0.15                         # gene letter sits just above arrow
    gap_y = -(band_h + 0.20)                  # gap (bp) sits BELOW the operon band
    badge_y = tag_y + 0.30                    # operon-id caption sits above letters
    # ONE row's data allocation. Constant (NOT ht-dependent) and equal to
    # _TRACK_YRANGE so the axis data-range below is exactly nrows*row_pitch: that
    # makes inches-per-data-unit identical across pages -> arrows are the same
    # physical thickness regardless of how many rows an operon wraps onto.
    row_pitch = _TRACK_YRANGE

    def rc(k):                                # (row, local column, baseline y) of gene k
        r = k // per
        return r, k - r * per, -r * row_pitch

    # operon grouping bands over maximal same-operon runs (neutral tint), drawn
    # first so the backbone and arrows sit on top; a run never crosses a row.
    i = 0
    while i < n:
        oid = members[i].get("operon_id") or ""
        r_i, c_i, y_i = rc(i)
        j = i
        while j + 1 < n and (members[j + 1].get("operon_id") or "") == oid \
                and (j + 1) // per == r_i:
            j += 1
        c_j = j - r_i * per
        if is_operon(oid):
            ax.add_patch(FancyBboxPatch(
                (c_i + off - 0.5, y_i - band_h), (c_j - c_i) + 1.0, 2 * band_h,
                boxstyle="round,pad=0.004,rounding_size=0.06",
                linewidth=1.0, edgecolor="#c3ccd6", facecolor="#eef1f6",
                alpha=0.95, zorder=0, mutation_aspect=0.5))
            if badge_per_operon and oid in badge_per_operon and r_i == 0:
                ax.text((c_i + c_j) / 2 + off, y_i + badge_y, badge_per_operon[oid],
                        ha="center", va="bottom", fontsize=8.6, fontstyle="italic",
                        color="#6b7280", zorder=4)
        i = j + 1

    # genome backbone behind the arrows -- one segment per row
    for r in range(nrows):
        y_r = -r * row_pitch
        if nrows > 1:
            last_col = (min((r + 1) * per, n) - 1) - r * per
            # extend 0.60 past the last gene, matching the 0.60 before the first
            # (symmetric backbone, same as the single-row case) so a full wrapped
            # row spans the same fraction left & right and the page stays symmetric
            bx0, bx1 = off - 0.60, last_col + off + 0.60
        else:
            bx0, bx1 = x0, x1
        ax.plot([bx0, bx1], [y_r, y_r], color="#b7bdc6", lw=1.6, zorder=1,
                solid_capstyle="round")

    # (wrapped rows read as one operon from the continuous gene numbering + the
    # "(N rows)" note in the title; no connector line needed.)

    # intergenic gaps (bp) under the backbone (within a row only)
    if show_gaps:
        for i in range(n - 1):
            r_i, c_i, y_i = rc(i)
            if (i + 1) // per != r_i:         # neighbour wraps to the next row
                continue
            a, b = members[i], members[i + 1]
            try:
                gap = int(b.get("start")) - int(a.get("end")) - 1
            except (TypeError, ValueError):
                continue
            ax.text(c_i + off + 0.5, y_i + gap_y, f"{gap:,} bp", ha="center",
                    va="top", fontsize=gap_fs, color="#7a7f87", zorder=2)

    entries = []
    for i, m in enumerate(members):
        r_i, c_i, y_i = rc(i)
        ax.add_patch(_gene_arrow(c_i + off, y_i, m.get("strand", "+"), gene_col[i], geom))
        tag = tag_letter(letter_offset + i)
        ax.text(c_i + off, y_i + tag_y, tag, ha="center", va="bottom",
                fontsize=tag_fs, fontweight="bold", color=gene_col[i], zorder=4)
        full = member_descriptor(m)             # source-prefixed, cleaned descriptor
        score = {"c2": m.get("c2"), "c3": m.get("c3"), "prelim": m.get("preliminary"),
                 "operon": m.get("operon_adjusted"), "final": m.get("confidence"),
                 "c3_hybrid": m.get("c3_hybrid"), "final_hybrid": m.get("operon_adjusted_hybrid"),
                 "feature_type": m.get("feature_type"),
                 "review": m.get("needs_review"),
                 "c1": m.get("c1"), "c4": m.get("c4"),
                 "review_reason": m.get("review_reason")}
        loc = ""
        try:
            loc = f"{int(m['start']):,}–{int(m['end']):,} ({m.get('strand', '+')})"
        except (TypeError, ValueError, KeyError):
            loc = ""
        entries.append((tag, full, score, loc, gene_col[i]))

    ax.set_xlim(x0 - 0.10, x1 + 0.10)
    top_extra = 0.42 if badge_per_operon else 0.0
    # total data range = EXACTLY nrows*row_pitch (see row_pitch note) so the
    # data->inch scale matches _arrow_geometry's nominal and arrows stay constant.
    y_top = ht + 0.48 + top_extra
    ax.set_ylim(y_top - nrows * row_pitch, y_top)
    # stash what pin_track_scale() needs to re-pin this track exactly after layout
    ax._track_scale = {"y_top": y_top, "nrows": nrows}
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ("top", "right", "bottom", "left"):
        ax.spines[s].set_visible(False)
    return entries


def _fmt_num(v):
    return None if v is None or (isinstance(v, float) and np.isnan(v)) else float(v)


def _score_breakdown(score) -> str:
    """'C3 0.65 · confidence = (prelim 0.78, operon +0.05, final 0.83)' from a
    score dict, omitting parts that are missing. A plain float -> '(0.83)'."""
    if score is None:
        return ""
    if isinstance(score, (int, float)):
        v = _fmt_num(score)
        return "" if v is None else f"({v:.2f})"
    c3 = _fmt_num(score.get("c3"))
    p = _fmt_num(score.get("prelim"))
    op = _fmt_num(score.get("operon"))
    f = _fmt_num(score.get("final"))
    parts = []
    if c3 is not None:
        parts.append(f"C3 {c3:.2f}")
    conf_bits = []
    if p is not None:
        conf_bits.append(f"prelim {p:.2f}")
    if p is not None and op is not None:
        conf_bits.append(f"operon {op - p:+.2f}")
    if f is not None:
        conf_bits.append(f"final {f:.2f}")
    if conf_bits:
        parts.append("confidence = (" + ", ".join(conf_bits) + ")")
    return "   ".join(parts)


def render_index(ax, entries, ncols: int = 1, fontsize: float = 10.5,
                 header: str | None = None, two_line: bool = False) -> None:
    """Render the gene-tag index into an axis. Full descriptors, never truncated.
    `entries` = [(tag, descriptor, score)], score a dict {c3,prelim,operon,final}
    or a float or None. two_line=True puts the score breakdown on an indented
    second line under each descriptor (used for the operon galleries)."""
    ax.axis("off")
    if not entries:
        return
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    top = 0.98
    if header:
        ax.text(0.0, 1.0, header, ha="left", va="top", fontsize=fontsize + 0.5,
                fontweight="bold", color=INK)
        top = 0.88
    lines_per = 2 if two_line else 1
    per_col = max(int(np.ceil(len(entries) / ncols)), 1)
    col_w = 1.0 / ncols
    line_h = min(top / (per_col * lines_per + 0.6), 0.16 if two_line else 0.30)
    tag_gap = 0.035 / ncols
    for k, entry in enumerate(entries):
        tag, desc, score = entry[0], entry[1], entry[2]
        c = k // per_col
        r = k % per_col
        x = c * col_w
        y = top - r * lines_per * line_h
        breakdown = _score_breakdown(score)
        if two_line:
            ax.text(x, y, tag, ha="left", va="top", fontsize=fontsize,
                    fontweight="bold", color=INK)
            ax.text(x + tag_gap, y, desc, ha="left", va="top", fontsize=fontsize,
                    color=INK)
            if breakdown:
                ax.text(x + tag_gap, y - line_h, breakdown, ha="left", va="top",
                        fontsize=fontsize - 1.2, color="#555555")
        else:
            suffix = f"  {breakdown}" if breakdown else ""
            ax.text(x, y, tag, ha="left", va="top", fontsize=fontsize,
                    fontweight="bold", color=INK)
            ax.text(x + tag_gap, y, f"{desc}{suffix}", ha="left", va="top",
                    fontsize=fontsize, color="#333333")


_TABLE_DESC_WRAP = 46   # descriptor wrap width (chars); long names wrap, never truncate
_TABLE_HEADER_UNITS = 2.4   # header takes 2 lines now (name + factor sub-line) + pad


def wrap_desc(desc: str, width: int = _TABLE_DESC_WRAP) -> list[str]:
    return textwrap.wrap(desc or "—", width=width) or ["—"]


def table_line_units(entries, desc_wrap: int = _TABLE_DESC_WRAP) -> float:
    """Total text lines a gene table needs (2-line header + wrapped descriptor
    lines), so a caller can size the table axis before drawing."""
    return _TABLE_HEADER_UNITS + sum(len(wrap_desc(e[1], desc_wrap)) for e in entries)


def _conflict_code(score) -> str:
    """Short 'conflict type' tag for the atlas table, read off the scored signals:
    EC  (C4<1 -> independent EC sources disagree),
    ambig (operon inference ambiguous -- mostly-uncharacterised operon context),
    desc (descriptor-consensus conflict). Em-dash when none detected."""
    codes = []
    c4 = _fmt_num(score.get("c4"))
    if c4 is not None and c4 < 0.999:
        codes.append("EC")
    reason = (score.get("review_reason") or "").lower()
    if "ambiguous" in reason:
        codes.append("ambig")
    if "descriptor" in reason and "conflict" in reason:
        codes.append("desc")
    return "+".join(codes) if codes else "—"


def _review_reason_short(score) -> str:
    """Readable short 'why' for the review flag, shown in the atlas 'review reason'
    column: the conflict / ambiguity / low-confidence triggers, condensed to a
    single line. Em-dash when the gene is not flagged."""
    parts = []
    reason = (score.get("review_reason") or "").lower()
    c4 = _fmt_num(score.get("c4"))
    if (c4 is not None and c4 < 0.999) or "ec conflict" in reason:
        parts.append("EC conflict")
    if "ambiguous" in reason:
        parts.append("ambiguous operon")
    if "descriptor" in reason and "conflict" in reason:
        parts.append("descriptor conflict")
    if "weak operon" in reason:
        parts.append("weak operon prob")
    if "low confidence" in reason:
        parts.append("low confidence")
    s = "; ".join(parts) if parts else "—"
    return s if len(s) <= 48 else s[:47] + "…"


def breakdown_col_layout(entries, fontsize: float, page_width_in: float,
                         desc_wrap: int, left: float = 0.0, right: float = 0.985) -> dict:
    """Compute the full-breakdown table's column x-positions (axis fractions) from
    the CONTENT of `entries`, laid out between the [left, right] axis-fraction
    bounds (default full axis; pass the arrow-backbone extent so the table sits
    within the SAME margins as the operon map). The fixed-content columns (gene,
    location, C1-C4, prelim, boost, final, review?) take only the width their
    values need; the two variable TEXT columns -- descriptor (left) and review
    reason (right) -- share all the leftover width so, left-aligned, they stretch
    the table across the bounds. Call ONCE per page over EVERY operon's entries and
    pass the result to each render_gene_table(col_layout=...) so all tables align."""
    fs, fsn, fsf = fontsize, fontsize - 0.4, fontsize - 2.1
    Win = page_width_in or 18.0
    def _wf(s, pt):
        return (len(str(s)) * 0.60 * pt / 72.0) / Win
    GAP = 0.24 / Win
    def _numw(h, sub):
        return max(_wf(h, fs), _wf(sub, fsf), _wf("+0.00", fsn))
    def _numw2(h, sub):        # holds an "adj/hyb" pair such as 0.00/0.00
        return max(_wf(h, fs), _wf(sub, fsf), _wf("0.00/0.00", fsn))
    wrapped = [wrap_desc(e[1], desc_wrap) for e in entries]
    tag_w = max([_wf("gene", fs)] + [_wf(e[0], fs) for e in entries])
    loc_w = max([_wf("location (bp)", fs)]
                + [_wf(e[3] if len(e) > 3 else "", fsn) for e in entries])
    type_w = max(_wf("type", fs), _wf("CDS/RNA", fsf), _wf("prophage", fsn))
    c1w, c2w = _numw("C1", "tool cov"), _numw("C2", "operon")
    c3w, c4w = _numw2("C3", "cons adj/hyb"), _numw("C4", "EC agree")
    prew = _numw("prelim", "C1×C4")
    opw = max(_numw("operon boost", "C2×C3 a/h"), _wf("+0.00/+0.00", fsn))
    finw = _numw2("final", "clip adj/hyb")
    revw = max(_wf("review?", fs), _wf("needs_review", fsf), _wf("yes", fsn))
    desc_nat = max([_wf("best_consensus_product_descriptor", fs)]
                   + [_wf(wl, fs) for w in wrapped for wl in w])
    reasons = [_review_reason_short(e[2]) if isinstance(e[2], dict) else "" for e in entries]
    reason_nat = max([_wf("review reason", fs)] + [_wf(r, fsn) for r in reasons])
    fixed = (tag_w + loc_w + type_w + c1w + c2w + c3w + c4w + prew + opw + finw + revw
             + 12 * GAP)
    avail, need = (right - left) - fixed, desc_nat + reason_nat
    if avail >= need:                                  # spare room -> share 60/40
        extra = avail - need
        w_desc, w_reason = desc_nat + 0.60 * extra, reason_nat + 0.40 * extra
    elif need > 0:                                      # too tight -> shrink both
        w_desc, w_reason = desc_nat * avail / need, reason_nat * avail / need
    else:
        w_desc = w_reason = 0.0
    x, pos = left, {}
    pos["x_tag"] = x; x += tag_w + GAP
    pos["x_desc"] = x; x += w_desc + GAP
    pos["x_loc"] = x; x += loc_w + GAP
    pos["x_type"] = x; x += type_w + GAP
    pos["x_c1"] = x; x += c1w + GAP
    pos["x_c2"] = x; x += c2w + GAP
    pos["x_c3"] = x; x += c3w + GAP
    pos["x_c4"] = x; x += c4w + GAP
    pos["x_pre"] = x; x += prew + GAP
    pos["x_op"] = x; x += opw + GAP
    pos["x_fin"] = x; x += finw + GAP
    pos["x_rev"] = x; x += revw + GAP
    pos["x_reason"] = x
    return pos


def _fit_desc_wrap(members_per, fontsize: float, page_width_in: float,
                   desc_wrap: int, left: float, right: float) -> int:
    """Largest descriptor wrap width (chars, <= desc_wrap) at which the descriptor
    column fits between the gene tag and the location column without spilling into
    it. Uses the same width model as breakdown_col_layout, applied BEFORE the table
    height and column layout are computed, so the wrap, the row count, and the
    columns all agree. Returns desc_wrap unchanged when there is already room."""
    fs, fsn, fsf = fontsize, fontsize - 0.4, fontsize - 2.1
    Win = page_width_in or 18.0
    def _wf(s, pt):
        return (len(str(s)) * 0.60 * pt / 72.0) / Win
    GAP = 0.24 / Win
    def _numw(h, sub):
        return max(_wf(h, fs), _wf(sub, fsf), _wf("+0.00", fsn))
    def _numw2(h, sub):
        return max(_wf(h, fs), _wf(sub, fsf), _wf("0.00/0.00", fsn))
    members = [m for ms in members_per for m in ms]
    if not members:
        return desc_wrap
    def _loc(m):
        try:
            return f"{int(m['start']):,}–{int(m['end']):,} ({m.get('strand', '+')})"
        except (TypeError, ValueError, KeyError):
            return ""
    tag_w = max([_wf("gene", fs)]
                + [_wf(str(i + 1), fs) for ms in members_per for i in range(len(ms))])
    loc_w = max([_wf("location (bp)", fs)] + [_wf(_loc(m), fsn) for m in members])
    type_w = max(_wf("type", fs), _wf("CDS/RNA", fsf), _wf("prophage", fsn))
    c1w, c2w = _numw("C1", "tool cov"), _numw("C2", "operon")
    c3w, c4w = _numw2("C3", "cons adj/hyb"), _numw("C4", "EC agree")
    prew = _numw("prelim", "C1×C4")
    opw = max(_numw("operon boost", "C2×C3 a/h"), _wf("+0.00/+0.00", fsn))
    finw = _numw2("final", "clip adj/hyb")
    revw = max(_wf("review?", fs), _wf("needs_review", fsf), _wf("yes", fsn))
    reasons = [_review_reason_short(m) for m in members]
    reason_nat = max([_wf("review reason", fs)] + [_wf(r, fsn) for r in reasons])
    fixed = (tag_w + loc_w + type_w + c1w + c2w + c3w + c4w + prew + opw + finw + revw
             + 12 * GAP)
    # width the descriptor column can take once the fixed columns and the reason
    # column (kept at its natural width) are placed
    w_desc = (right - left) - fixed - reason_nat
    char_w = 0.60 * fs / 72.0 / Win
    n = int(w_desc / char_w) if char_w > 0 else desc_wrap
    return max(20, min(desc_wrap, n))


def render_gene_table(ax, entries, fontsize: float = 7.6, show_location: bool = True,
                      show_scores: bool = True, desc_wrap: int = _TABLE_DESC_WRAP,
                      full_breakdown: bool = False, page_width_in: float = None,
                      col_layout: dict = None) -> None:
    """Draw a compact gene table in `ax`. Columns: gene tag (coloured to match
    its arrow), full descriptor WRAPPED across lines (never truncated), genomic
    location, the score breakdown, and a review flag. Each score column carries
    a small sub-header naming what it is a function of:
      C3 (operon coherence) · prelim = C1·C4 · operon boost = f(C2,C3) ·
      final = min(1, prelim+context) · review? (yes only when the operon
      context LOWERED the score by >= 0.1; smaller drops and positive/zero
      boosts are never flagged).
    `entries` are the 5-tuples from draw_gene_track: (tag, descriptor, score,
    location, colour); score is a dict with c3/prelim/operon/final/review."""
    ax.axis("off")
    if not entries:
        return
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    have_scores = show_scores and any(
        isinstance(e[2], dict) and _fmt_num(e[2].get("final")) is not None for e in entries)
    have_loc = show_location and any(len(e) > 3 and e[3] for e in entries)
    # review flag is driven ONLY by the operon context LOWERING the score by a
    # material amount (>= 0.1): smaller drops are not significant, and a
    # positive/zero operon boost never triggers review.
    have_review = show_scores and any(
        isinstance(e[2], dict) and _fmt_num(e[2].get("prelim")) is not None
        and _fmt_num(e[2].get("operon")) is not None for e in entries)
    fs, fsn, fsf = fontsize, fontsize - 0.4, fontsize - 2.1
    wrapped = [wrap_desc(e[1], desc_wrap) for e in entries]

    # column x-positions.
    if full_breakdown:
        # DYNAMIC, all-LEFT-aligned layout (see breakdown_col_layout). A page-level
        # col_layout is passed in so every stacked operon table on the page shares
        # the SAME columns; fall back to a per-table layout if none was given.
        pos = col_layout or breakdown_col_layout(entries, fontsize, page_width_in, desc_wrap)
        x_tag, x_desc, x_loc = pos["x_tag"], pos["x_desc"], pos["x_loc"]
        x_type = pos["x_type"]
        x_c1, x_c2, x_c3, x_c4 = pos["x_c1"], pos["x_c2"], pos["x_c3"], pos["x_c4"]
        x_pre, x_op, x_fin = pos["x_pre"], pos["x_op"], pos["x_fin"]
        x_rev, x_reason, x_conf = pos["x_rev"], pos["x_reason"], None
    else:
        x_tag, x_desc = 0.0, 0.03
        x_c1 = x_c4 = x_conf = x_reason = x_type = None
        x_loc, x_c2, x_c3, x_pre, x_op, x_fin, x_rev = 0.335, 0.485, 0.565, 0.652, 0.738, 0.838, 0.935

    units = _TABLE_HEADER_UNITS + sum(len(w) for w in wrapped)
    line_h = 0.992 / units
    y = 0.992

    def cell(x, text, bold=False, color="#222222", size=fs, yy=None, ha="left"):
        ax.text(x, y if yy is None else yy, text, ha=ha, va="top", fontsize=size,
                fontweight="bold" if bold else "normal", color=color)

    cell(x_tag, "gene", bold=True)
    cell(x_desc, "best_consensus_product_descriptor", bold=True)
    if have_loc:
        cell(x_loc, "location (bp)", bold=True)
    if full_breakdown and x_type is not None:
        cell(x_type, "type", bold=True)
    if have_scores:
        if full_breakdown:
            cell(x_c1, "C1", bold=True)
        cell(x_c2, "C2", bold=True); cell(x_c3, "C3", bold=True)
        if full_breakdown:
            cell(x_c4, "C4", bold=True)
        cell(x_pre, "prelim", bold=True)
        cell(x_op, "operon boost", bold=True); cell(x_fin, "final", bold=True)
    if have_review:
        cell(x_rev, "review?", bold=True)
        if full_breakdown:
            cell(x_reason, "review reason", bold=True)
    # factor sub-header line (smaller, grey) -- names each column's meaning so the
    # operon-boost trail is self-explanatory (C2 = operon probability GATE, C3 =
    # cross-genome conservation; boost = C2·C3, penalty = C2·conflict)
    ysub = y - line_h
    if full_breakdown:
        subs = [(x_type, "feature"), (x_c1, "tool cov"), (x_c2, "operon"), (x_c3, "cons adj/hyb"),
                (x_c4, "EC agree"), (x_pre, "C1×C4"), (x_op, "C2×C3 a/h"), (x_fin, "clip adj/hyb")]
    else:
        subs = [(x_c2, "operon prob."), (x_c3, "conservation"), (x_pre, "C1 × C4"),
                (x_op, "C2×C3 (boost)"), (x_fin, "min(1, prelim+ctx)")]
    for x, sub in subs:
        if have_scores and x is not None:
            cell(x, sub, size=fsf, color=_DETAIL_PURPLE, yy=ysub)
    if have_review:
        cell(x_rev, ("needs_review" if full_breakdown else "operon lowered ≥0.1"),
             size=fsf, color=_DETAIL_PURPLE, yy=ysub)
        if full_breakdown:
            cell(x_reason, "why flagged", size=fsf, color=_DETAIL_PURPLE, yy=ysub)
    y -= line_h * _TABLE_HEADER_UNITS

    for entry, wlines in zip(entries, wrapped):
        tag, desc, score = entry[0], entry[1], entry[2]
        loc = entry[3] if len(entry) > 3 else ""
        color = entry[4] if len(entry) > 4 else INK
        y0 = y
        ax.text(x_tag, y0, tag, ha="left", va="top", fontsize=fs, fontweight="bold",
                color=color)
        for li, wl in enumerate(wlines):
            ax.text(x_desc, y0 - li * line_h, wl, ha="left", va="top", fontsize=fs,
                    color=INK)
        if have_loc and loc:
            ax.text(x_loc, y0, loc, ha="left", va="top", fontsize=fsn, color="#555555")
        if full_breakdown and x_type is not None and isinstance(score, dict):
            fty = str(score.get("feature_type") or "").strip()
            if fty:
                # non-CDS features (rna, prophage) flagged in a warm tone so the
                # reader can tell at a glance a row is not a protein-coding gene.
                tcol = "#555555" if fty.lower() == "cds" else "#a86400"
                ax.text(x_type, y0, fty, ha="left", va="top", fontsize=fsn, color=tcol)
        if have_scores and isinstance(score, dict):
            c2 = _fmt_num(score.get("c2")); c3 = _fmt_num(score.get("c3"))
            p = _fmt_num(score.get("prelim"))
            op = _fmt_num(score.get("operon")); f = _fmt_num(score.get("final"))
            c1 = _fmt_num(score.get("c1")); c4 = _fmt_num(score.get("c4"))
            # colour a driver RED when it sits below the 0.5 neutral line -- that is
            # exactly what pulls the operon boost negative, so the cause is visible.
            if full_breakdown and c1 is not None:
                ax.text(x_c1, y0, f"{c1:.2f}", va="top", fontsize=fsn,
                        color=(LOWERED if c1 < 0.5 else "#333333"))
            if c2 is not None:
                ax.text(x_c2, y0, f"{c2:.2f}", va="top", fontsize=fsn,
                        color=(LOWERED if c2 < 0.5 else "#333333"))
            if c3 is not None:
                c3h = _fmt_num(score.get("c3_hybrid"))
                txt = f"{c3:.2f}/{c3h:.2f}" if (full_breakdown and c3h is not None) else f"{c3:.2f}"
                ax.text(x_c3, y0, txt, va="top", fontsize=fsn,
                        color=(LOWERED if c3 < 0.5 else "#333333"))
            if full_breakdown and c4 is not None:
                # C4 = EC-conflict clearance; < 1 means an EC conflict was applied.
                ax.text(x_c4, y0, f"{c4:.2f}", va="top", fontsize=fsn,
                        color=(LOWERED if c4 < 0.999 else "#333333"))
            if p is not None:
                ax.text(x_pre, y0, f"{p:.2f}", va="top", fontsize=fsn, color="#333333")
            if p is not None and op is not None:
                if full_breakdown:
                    # RAW operon context term = C2·C3 (boost-only, c2-gated), shown
                    # adj/hyb -- the actual model term the header names, NOT
                    # final-prelim (which would just echo the final column, and hides
                    # any headroom lost to the clip at 1).
                    c2v = _fmt_num(score.get("c2"))
                    c3a = _fmt_num(score.get("c3"))
                    c3h = _fmt_num(score.get("c3_hybrid"))
                    ba = (max(0.0, c2v) * max(0.0, c3a)) if (
                        c2v is not None and c3a is not None) else None
                    bh = (max(0.0, c2v) * max(0.0, c3h)) if (
                        c2v is not None and c3h is not None) else None
                    if ba is not None and bh is not None:
                        txt, bcol = f"{ba:+.2f}/{bh:+.2f}", (RAISED if max(ba, bh) > 0 else "#333333")
                    elif ba is not None:
                        txt, bcol = f"{ba:+.2f}", (RAISED if ba > 0 else "#333333")
                    else:
                        txt = None
                    if txt:
                        ax.text(x_op, y0, txt, va="top", fontsize=fsn, color=bcol)
                else:
                    d = op - p                        # compact report: effective delta
                    ax.text(x_op, y0, f"{d:+.2f}", va="top", fontsize=fsn,
                            color=(RAISED if d > 0 else (LOWERED if d < 0 else "#333333")))
            if f is not None:
                fh = _fmt_num(score.get("final_hybrid"))
                txt = f"{f:.2f}/{fh:.2f}" if (full_breakdown and fh is not None) else f"{f:.2f}"
                ax.text(x_fin, y0, txt, va="top", fontsize=fsn,
                        fontweight="bold", color="#111111")
        if have_review and isinstance(score, dict):
            if full_breakdown:
                # the ACTUAL scored review flag (EC conflict / operon ambiguity /
                # low confidence), not the operon-lowered heuristic -- which never
                # fires under boost-only -- plus a readable reason column.
                nr = str(score.get("review") or "").strip().lower()
                flagged = nr in ("yes", "true", "1")
                ax.text(x_rev, y0, "yes" if flagged else "no", va="top", fontsize=fsn,
                        fontweight="bold" if flagged else "normal",
                        color=(LOWERED if flagged else "#5a9e6f"))
                ax.text(x_reason, y0, _review_reason_short(score), va="top",
                        fontsize=fsn, color=(LOWERED if flagged else "#8a8f97"))
            else:
                p_r = _fmt_num(score.get("prelim")); op_r = _fmt_num(score.get("operon"))
                if p_r is not None and op_r is not None:
                    if (op_r - p_r) <= -0.1:  # operon context materially LOWERED the score (>= 0.1)
                        ax.text(x_rev, y0, "yes", va="top", fontsize=fsn, fontweight="bold",
                                color=LOWERED)
                    else:
                        ax.text(x_rev, y0, "no", va="top", fontsize=fsn, color="#5a9e6f")
        y -= len(wlines) * line_h


def render_operon_page(outpath, blocks, *, org_label, suptitle, run_root=None,
                       provenance=None, note=None, footer_sources=None,
                       fig_width=22.0, table_fs=13.0, desc_wrap=60,
                       per_row=None, min_span=None, dpi=170, badges=None,
                       legend_handles=None) -> None:
    """Render ONE page of operon blocks (arrow track + full C1-C4 breakdown table)
    in the UNIFIED format shared by the full-genome atlas and the representative
    report galleries, so their images look identical:
      * deterministic absolute placement -- every arrow row is _TRACK_ROW_IN tall,
        so arrows are a constant physical size no matter how an operon wraps;
      * ONE page-level column layout (breakdown_col_layout over EVERY block's
        entries) so all stacked tables align their columns;
      * a CENTERED content column (x0 = (1-xw)/2) -> equal left/right page margins
        and the table sits within the same margins as the arrow map;
      * a font hierarchy scaled off `table_fs` (one knob for the whole page).
    `blocks` = list of {"members": [...], "title": "..."}.  Returns nothing; writes
    the PNG to `outpath`."""
    n = len(blocks)
    if n == 0:
        return
    members_per = [b["members"] for b in blocks]

    # Page geometry, computed up front so the descriptor wrap can be fit to the
    # real column width -- this keeps the wrap, the table row count, and the column
    # layout in agreement so long descriptors never spill into the location column.
    xw = 0.945
    x0 = (1.0 - xw) / 2.0                  # CENTERED axis -> equal left/right margins
    page_width_in = fig_width * xw
    span_ref = min_span or max((len(m) for m in members_per), default=1)
    left_frac = 0.10 / (span_ref + 0.40)
    desc_wrap = _fit_desc_wrap(members_per, table_fs, page_width_in, desc_wrap,
                               left_frac, 1.0 - left_frac)

    def _nrows(m):
        return (-(-len(m) // per_row)) if (per_row and len(m) > per_row) else 1
    nrows_per = [_nrows(m) for m in members_per]
    units_per = [member_table_units(m, desc_wrap) for m in members_per]
    track_h = [_TRACK_ROW_IN * nr for nr in nrows_per]
    line_in = round(table_fs / 72.0 * 1.42, 3)
    table_h = [max(units_per[i], 1.0) * line_in for i in range(n)]

    # one knob (table_fs) drives the whole font hierarchy so it scales cleanly
    title_fs = table_fs + 1.0          # per-operon "N-gene operon | ..." title
    tag_fs = table_fs + 1.0            # gene numbers above the arrows
    gap_fs = table_fs - 1.5            # intergenic "N bp" labels
    prov_fs = table_fs - 0.5           # provenance line
    note_fs = table_fs - 1.0           # method note / legend
    foot_fs = table_fs - 2.0           # sources footer
    sup_fs = table_fs + 4.0            # page suptitle
    org_fs = table_fs + 2.0            # organism scientific name

    has_legend = bool(legend_handles)
    _HEADER_IN = 2.30 if has_legend else 1.90
    _TITLE_GAP, _INTRA_GAP, _INTER_GAP, _FOOTER_IN = 0.46, 0.12, 0.38, 0.85
    block_in = [_TITLE_GAP + track_h[i] + _INTRA_GAP + table_h[i] for i in range(n)]
    H = _HEADER_IN + sum(block_in) + _INTER_GAP * (n - 1) + _FOOTER_IN

    # The arrow backbone spans axis x = -0.60..span-0.40 over an xlim of span+0.40,
    # i.e. axis fractions [left_frac, 1-left_frac]. Lay the title, the table and the
    # arrows all between those SAME mirrored bounds so the ink is symmetric within
    # the axis -> the page crops to EQUAL left/right margins and the table sits
    # within the same margins as the operon map.
    fig = plt.figure(figsize=(fig_width, H))
    track_axes, tables = [], []
    cur_top = H - _HEADER_IN
    for i, b in enumerate(blocks):
        members = members_per[i]
        t_top = cur_top - _TITLE_GAP
        t_bot = t_top - track_h[i]
        axT = fig.add_axes([x0, t_bot / H, xw, track_h[i] / H])
        tab_top = t_bot - _INTRA_GAP
        tab_bot = tab_top - table_h[i]
        axTab = fig.add_axes([x0, tab_bot / H, xw, table_h[i] / H])
        # hide the axes' background rectangles so savefig's tight crop bounds to the
        # actual INK, not the full-width axis patch
        axT.patch.set_visible(False)
        axTab.patch.set_visible(False)
        track_axes.append(axT)
        wrap = per_row if (per_row and len(members) > per_row) else None
        entries = draw_gene_track(axT, members, show_gaps=True, min_span=min_span,
                                  per_row=wrap, badge_per_operon=badges,
                                  tag_fs=tag_fs, gap_fs=gap_fs)
        # inset the title to left_frac (not the axis edge) so it doesn't anchor the
        # far-left margin past where the arrow map / table start
        ttl = axT.set_title(b["title"], fontsize=title_fs, pad=6, loc="left",
                            fontweight="bold")
        ttl.set_x(left_frac)
        tables.append((axTab, entries))
        cur_top = tab_bot - _INTER_GAP

    # single page-level column layout from EVERY block's entries -> aligned tables
    all_entries = [e for _, es in tables for e in es]
    col_layout = breakdown_col_layout(all_entries, table_fs, page_width_in, desc_wrap,
                                      left=left_frac, right=1.0 - left_frac)
    for axTab, entries in tables:
        render_gene_table(axTab, entries, fontsize=table_fs, desc_wrap=desc_wrap,
                          full_breakdown=True, col_layout=col_layout)
    for axT in track_axes:
        pin_track_scale(axT)

    fh = H
    st = fig.suptitle(suptitle, y=1 - 0.48 / fh, fontweight="bold", fontsize=sup_fs)
    extras = [st]
    draw_organism_line(fig, org_label, y=1 - 1.02 / fh, fontsize=org_fs)
    if provenance is not None:
        set_provenance(provenance)
    pv = draw_provenance_line(fig, 1 - 1.46 / fh, fontsize=prov_fs)
    if pv is not None:
        extras.append(pv)
    if has_legend:
        lg = fig.legend(handles=legend_handles, loc="upper center",
                        ncol=len(legend_handles), bbox_to_anchor=(0.5, 1 - 1.86 / fh),
                        fontsize=note_fs, frameon=False, handlelength=1.4)
        extras.append(lg)
    fig._report_extra_artists = getattr(fig, "_report_extra_artists", []) + extras
    if note:
        draw_method_note(fig, note, fontsize=note_fs)
    if footer_sources and run_root is not None:
        draw_sources_footer(fig, run_root, footer_sources, fontsize=foot_fs)
    savefig(fig, Path(outpath), dpi=dpi)


def operon_members_informative(members_in_order: str, min_informative: int | None = None):
    """Judge an operon (its "a -> b -> ..." label string) by member informativeness.
    Returns True when the operon has >=2 members and either ALL are named-function
    (min_informative=None) or at least `min_informative` are. Uses the same
    is_uninformative gate as scoring, so "hypothetical protein -> hypothetical
    protein" and blank "-" operons are treated as not-named-function."""
    members = [m.strip() for m in (members_in_order or "").split(" -> ") if m.strip()]
    if len(members) < 2:
        return False
    informative = [(m not in ("-", "")) and not is_uninformative(m) for m in members]
    if min_informative is None:
        return all(informative)
    return sum(informative) >= min_informative


def operon_to_members(op_row) -> list[dict]:
    """Convert a build_operons() row into draw_gene_track members, carrying the
    per-gene score breakdown (C3, preliminary, operon-adjusted, final)."""
    out = []
    labels = op_row["member_labels"]

    def col(name):
        return op_row[name] if name in op_row else [None] * len(labels)

    confs, prelims = col("confidences"), col("preliminaries")
    opadj, c3s, srcs = col("operon_adjusteds"), col("c3s"), col("sources")
    opadj_h, c3s_h = col("operon_adjusteds_hybrid"), col("c3s_hybrid")
    c2s = col("c2s")
    c1s, c4s = col("c1s"), col("c4s")
    nrs = col("needs_reviews")
    ftys = col("feature_types")
    reasons = col("review_reasons")
    for i in range(len(labels)):
        out.append({
            "start": op_row["starts"][i],
            "end": op_row["ends"][i],
            "strand": op_row["strands"][i],
            "label": labels[i],
            "source": srcs[i] if i < len(srcs) else None,
            "feature_type": ftys[i] if i < len(ftys) else None,
            "needs_review": nrs[i] if i < len(nrs) else None,
            "review_reason": reasons[i] if i < len(reasons) else None,
            "operon_id": op_row["operon_id"],
            "confidence": confs[i] if i < len(confs) else None,
            "preliminary": prelims[i] if i < len(prelims) else None,
            "operon_adjusted": opadj[i] if i < len(opadj) else None,
            "operon_adjusted_hybrid": opadj_h[i] if i < len(opadj_h) else None,
            "c3": c3s[i] if i < len(c3s) else None,
            "c3_hybrid": c3s_h[i] if i < len(c3s_h) else None,
            "c2": c2s[i] if i < len(c2s) else None,
            "c1": c1s[i] if i < len(c1s) else None,
            "c4": c4s[i] if i < len(c4s) else None,
        })
    return out
