#!/usr/bin/env python3
"""reportfig_lib.py - shared foundation for the post-scoring report figures.

Reads only finished per-organism scoring outputs and the depot operon reference,
never writes into scoring folders, and provides the style, layout, data loading,
gene-track and operon-page drawing used by every report figure (PNG plus TSV).
"""
from __future__ import annotations

import os
import re
import sys
import textwrap
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # headless / SLURM
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyBboxPatch, Polygon, Rectangle  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

# Uses scoring's c3_lib descriptor cleaning and uninformative gate so gene identity matches scoring.
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

# ---- palette ----
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
VERMILLION = RED   # aliases used by figure code
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

# Confidence tiers, best to worst, with fixed colours; non-coding genes are grey.
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

# Fixed colours for the confidence components.
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

# Serif preference: Times New Roman, then metric-compatible clones, then DejaVu Serif.
_SERIF = ["Times New Roman", "Nimbus Roman No9 L", "Nimbus Roman",
          "Liberation Serif", "Tinos", "Times", "DejaVu Serif"]

# Sans-serif preference: Helvetica/Arial or their clones, then DejaVu Sans.
_SANS = ["Helvetica", "Nimbus Sans", "Arial", "Liberation Sans",
         "TeX Gyre Heros", "DejaVu Sans"]

# Font files registered explicitly because matplotlib's default cache omits these directories.
_FONT_FILES = [
    "/usr/share/fonts/urw-base35/NimbusSans-Regular.otf",
    "/usr/share/fonts/urw-base35/NimbusSans-Bold.otf",
    "/usr/share/fonts/urw-base35/NimbusSans-Italic.otf",
    "/usr/share/fonts/urw-base35/NimbusSans-BoldItalic.otf",
    "/usr/share/fonts/liberation-sans/LiberationSans-Regular.ttf",
    "/usr/share/fonts/liberation-sans/LiberationSans-Bold.ttf",
    "/usr/share/fonts/liberation-sans/LiberationSans-Italic.ttf",
    # Times-metric serif clones.
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
    # mathtext italicises scientific names in a non-bold serif.
    "mathtext.fontset": "custom",
    "mathtext.rm": "serif",
    "mathtext.it": "serif:italic",
    "mathtext.bf": "serif:bold",
    "mathtext.default": "regular",
}


def _register_fonts() -> None:
    """Registers the bundled font files with matplotlib, skipping missing ones."""
    from matplotlib import font_manager as fm
    for p in _FONT_FILES:
        try:
            if Path(p).exists():
                fm.fontManager.addfont(p)
        except Exception:
            pass


def apply_style() -> None:
    """Registers fonts and applies the shared matplotlib rcParams."""
    _register_fonts()
    plt.rcParams.update(_STYLE)


# ---- layout helpers ----
def panel_letter(ax, letter: str) -> None:
    """Draws "(letter)" at the top-left corner, below and left of the panel title."""
    ax.text(-0.065, 1.012, f"({letter})", transform=ax.transAxes,
            fontsize=15, fontweight="bold", va="bottom", ha="left",
            color=INK, clip_on=False)


def set_title(ax, text: str) -> None:
    """Sets a bold centred title, padded to clear the panel letter."""
    ax.set_title(text, pad=26, loc="center", fontweight="bold")


# ---- provenance line ----
# Set once per run with set_provenance(); drawn under every figure title.
_PROVENANCE: str | None = None


def set_provenance(text: str | None) -> None:
    """Sets the provenance line drawn under every figure title this run."""
    global _PROVENANCE
    _PROVENANCE = text


def provenance_text(pool_organisms: int, stats=None, leave_one_out: bool = False) -> str:
    """Returns the provenance text describing the operon-database pool.

    stats is an aggregate_pool_stats() dict or a bare total-gene int; leave_one_out
    notes that the candidate genome was excluded from the pool."""
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
    """Draws the provenance line at figure-fraction y and returns the artist, or None if unset."""
    if not _PROVENANCE:
        return None
    return fig.text(0.5, y, _PROVENANCE, ha="center", va="top", fontsize=fontsize,
                    color=_DETAIL_PURPLE, style="italic")


def pool_total_genes(run_root, organisms) -> int:
    """Returns the total gene count over the organisms' scored output, skipping missing ones."""
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
    """Lays out the figure and reserves a title band of `band` inches.

    Draws the suptitle, an optional italic organism line and the provenance line."""
    if not suptitle and not organism:
        fig.tight_layout(h_pad=h_pad, w_pad=w_pad)
        return
    # tight_layout(rect=) ignores the band when a gene-track axis is present, so
    # the band is reserved afterwards with subplots_adjust.
    fh = max(fig.get_figheight(), 3.0)
    b = band if band is not None else (1.55 if organism else 1.4)
    if _PROVENANCE:
        b += 0.30  # headroom for the provenance line
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
    """Saves the figure cropped to its content plus any registered title/footer artists."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    extra = getattr(fig, "_report_extra_artists", None)
    if extra:
        # bbox_inches="tight" remaps figure-fraction artists, so the crop box is
        # computed once, including the suptitle that get_tightbbox would drop.
        r = fig.canvas.get_renderer()
        fig.draw(r)
        arts = list(extra)
        st = getattr(fig, "_suptitle", None)
        if st is not None and st not in arts:
            arts.append(st)
        # An explicit Bbox gets no margin, so it is padded here.
        bb = fig.get_tightbbox(r, bbox_extra_artists=arts).padded(0.40)
        fig.savefig(path, dpi=dpi, bbox_inches=bb, facecolor="white")
    else:
        fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"[reportfig] wrote {path.name}", file=sys.stderr)


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    """Writes df as a TSV, creating the parent directory."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)
    print(f"[reportfig] wrote {path.name}", file=sys.stderr)


def short_desc(desc: str, maxlen: int = 26) -> str:
    """Returns the cleaned descriptor truncated to maxlen with an ellipsis."""
    d = clean_descriptor(desc) or (desc or "")
    d = re.sub(r"\s+", " ", d).strip()
    return d if len(d) <= maxlen else d[: maxlen - 1] + "…"


def short_organism(org: str) -> str:
    """Returns a short name, e.g. 'Escherichia_coli_str._K-12...GCF_x' -> 'E. coli str. K-12'."""
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
    """Returns a mathtext organism label with the scientific name italic and strain upright.

    The assembly accession, if present, is appended in parentheses."""
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
    """Returns [(text, is_italic)] segments of an organism name for offsetbox rendering."""
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
    """Returns an HPacker of the organism segments with real italics."""
    from matplotlib.offsetbox import TextArea, HPacker
    segs = organism_segments(org)
    boxes = [TextArea(t, textprops=dict(
                fontstyle="italic" if it else "normal", fontweight="normal",
                fontsize=fontsize, color=INK)) for t, it in segs if t]
    return HPacker(children=boxes, align="baseline", pad=0, sep=fontsize * 0.30)


def draw_organism_line(fig, org: str, y: float = 0.945, fontsize: float = 15) -> None:
    """Draws the non-bold organism name centred at figure-fraction y as a Text artist."""
    t = fig.text(0.5, y, italic_organism(org), ha="center", va="top",
                 fontsize=fontsize, fontweight="normal", color=INK)
    fig._report_extra_artists = getattr(fig, "_report_extra_artists", []) + [t]


def draw_title(fig, description: str, org: str, fontsize_desc: float = 16,
               fontsize_org: float = 13.5, y: float = 0.998) -> None:
    """Draws a bold description over the organism line as one top-centred offsetbox."""
    from matplotlib.offsetbox import TextArea, VPacker, AnnotationBbox
    desc = TextArea(description, textprops=dict(fontweight="bold",
                    fontsize=fontsize_desc, color=INK))
    stack = VPacker(children=[desc, _organism_hpacker(org, fontsize_org)],
                    align="center", pad=0, sep=fontsize_org * 0.5)
    ab = AnnotationBbox(stack, (0.5, y), xycoords="figure fraction", frameon=False,
                        box_alignment=(0.5, 1.0))
    fig.add_artist(ab)
    fig._report_extra_artists = getattr(fig, "_report_extra_artists", []) + [ab]


# ---- data loading (finished scoring outputs only) ----
_CONF_FINAL = "scoring/scored-labeled-genes-confidence-final.tsv"
_LABELED = "labeling/labeled-genes.tsv"
_NOT_IN_OPERON = "NOT_IN_AN_OPERON"

# reorganize_outputs.py moves scoring/ and labeling/ into per-tool-phased-output/
# per organism, so a run can hold both layouts; files are resolved in either.
_PTP = "per-tool-phased-output"


def _resolve_organism_file(base, rel):
    """Returns the path of rel under an organism dir in either layout, or None."""
    base = Path(base)
    for cand in (base / rel, base / _PTP / rel):
        if cand.is_file():
            return cand
    return None


def _require_organism_file(base, rel):
    """Returns the path of rel in either layout, raising FileNotFoundError if absent."""
    p = _resolve_organism_file(base, rel)
    if p is None:
        raise FileNotFoundError(
            f"{rel} not found for {Path(base).name} in either layout "
            f"({rel} or {_PTP}/{rel})"
        )
    return p


def is_operon(oid) -> bool:
    """Returns True only for a real operon id (one starting with 'operon_')."""
    return bool(oid) and str(oid).startswith("operon_")

_NUMERIC_COLS = [
    "c1_score", "c2_score_from_operon_probability", "c3_score", "c4_score",
    "preliminary_confidence_c1_c4", "final_confidence_operon_context",
    "confidence_score", "operon_member_count", "operon_gene_position_in_operon",
]
_GENE_ID_RE = re.compile(r"^(.*)_(\d+)([+-])(\d+)$")


def parse_contig(gene_id: str) -> str:
    """Returns the contig part of a gene_id, or the id itself when it does not parse."""
    m = _GENE_ID_RE.match(gene_id or "")
    return m.group(1) if m else (gene_id or "")


def discover_organisms(run_root: Path) -> list[str]:
    """Returns organism stems under the run that have a confidence-final scoring file."""
    run_root = Path(run_root)
    out = set()
    # pre-reorganize:  <organism>/scoring/...
    for p in sorted(run_root.glob("*/" + _CONF_FINAL)):
        out.add(p.parents[1].name)
    # post-reorganize: <organism>/per-tool-phased-output/scoring/...
    for p in sorted(run_root.glob(f"*/{_PTP}/" + _CONF_FINAL)):
        out.add(p.parents[2].name)
    return sorted(out)


def _coerce_numeric(df: pd.DataFrame) -> pd.DataFrame:
    """Converts the known score columns to numbers in place."""
    for c in _NUMERIC_COLS:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def load_organism_genes(run_root: Path, organism: str) -> pd.DataFrame:
    """Returns one row per gene: confidence-final scores and operon info joined to coordinates.

    Adds clean_desc, uninformative, in_operon and contig."""
    run_root = Path(run_root)
    base = run_root / organism
    conf = pd.read_csv(_require_organism_file(base, _CONF_FINAL), sep="\t", dtype=str,
                       keep_default_na=False, engine="python")
    conf = _coerce_numeric(conf)

    coords_cols = ["feature_id", "gene_id", "gene_start", "gene_end", "RAST_strand"]
    lab = pd.read_csv(_require_organism_file(base, _LABELED), sep="\t", dtype=str,
                      keep_default_na=False, engine="python",
                      usecols=lambda c: c in coords_cols)
    for c in ("gene_start", "gene_end"):
        if c in lab.columns:
            lab[c] = pd.to_numeric(lab[c], errors="coerce")
    df = conf.merge(lab, on="feature_id", how="left")

    df["contig"] = df.get("gene_id", "").map(parse_contig)
    df["clean_desc"] = df["best_consensus_product_descriptor"].map(clean_descriptor)
    df["uninformative"] = df["best_consensus_product_descriptor"].map(is_uninformative)
    # Only ids starting with "operon_" are real operons.
    oid = df.get("operon_id", pd.Series([""] * len(df))).astype(str)
    df["in_operon"] = oid.str.startswith("operon_")
    df["organism"] = organism
    return df


def load_all_genes(run_root: Path, organisms: list[str] | None = None) -> pd.DataFrame:
    """Returns the concatenated gene tables of the organisms, skipping ones that fail to load."""
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


# ---- operons with ordered members and coordinates ----
def build_operons(genes: pd.DataFrame) -> pd.DataFrame:
    """Returns one row per operon with ordered members, coordinates and per-gene scores.

    members_in_order ("label1 -> label2 -> ...") is the join key of the depot operon database."""
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
            "tiers": g.get("confidence_tier", pd.Series([""] * len(g))).tolist(),
        })
    return pd.DataFrame(rows)


# ---- depot pangenome operon recurrence ----
def load_operon_recurrence(depot_db: Path, restrict_to=None) -> dict[str, dict]:
    """Returns {members_in_order: {label_frequency, organism_count, organisms}} from the depot DB.

    restrict_to, when non-empty, limits each operon's organisms to that set."""
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
        # Recurrence counts distinct organisms, so duplicates in the DB are dropped.
        orgs = sorted({x for x in (r.get("organisms", "") or "").split("|") if x})
        if keep is not None:
            orgs = [o for o in orgs if o in keep]
        n = len(orgs)
        rec = {
            # A scoped count cannot split the DB-wide frequency, so it uses the organism count.
            "label_frequency": (n if keep is not None
                                else int(r.get("fingerprint_label_frequency", 0) or 0)),
            "organism_count": n,
            "organisms": orgs,
        }
        # Keeps the record with most organisms when the same text appears under several hashes.
        prev = out.get(mio)
        if prev is None or rec["organism_count"] > prev["organism_count"]:
            out[mio] = rec
    return out


# Default depot locations of the operon recurrence DB and the OCC reference.
DEFAULT_OPERON_DB = Path(os.environ.get("MARGIE_OPERON_FINGERPRINT_DB", "operon-fingerprint-database-label-ordered.tsv"))
DEFAULT_OCC_REFERENCE = Path(os.environ.get("MARGIE_OCC_REFERENCE", "occ_reference.pkl"))


def load_occ_organisms(occ_reference=DEFAULT_OCC_REFERENCE):
    """Returns the set of organisms in the OCC reference (the pool C3 was scored against).

    Returns None when the reference is unreadable, so callers fall back to the run's organisms."""
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
    """Returns {organism: {stat: int}} from the OCC reference's .genome_stats.tsv sidecar.

    Returns an empty dict when the sidecar is absent or unreadable."""
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
    """Returns the pool stats summed over organisms, plus n_genomes (how many had stats)."""
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
    """Returns a footer path: run-relative tagged [output], else absolute tagged [host]."""
    p, run_root = Path(path), Path(run_root)
    try:
        return f"{p.relative_to(run_root)}  [output]"
    except ValueError:
        return f"{p}  [host]"


def organism_source_lines(run_root, organism: str, coords: bool = False,
                          operon_db=None) -> list[str]:
    """Returns the footer source list for a per-organism figure, using the layout actually read."""
    run_root = Path(run_root)
    _base = run_root / organism
    lines = [rel_or_host(_resolve_organism_file(_base, _CONF_FINAL)
                         or _base / _CONF_FINAL, run_root)]
    if coords:
        lines.append(rel_or_host(_resolve_organism_file(_base, _LABELED)
                                 or _base / _LABELED, run_root))
    if operon_db is not None:
        lines.append(rel_or_host(operon_db, run_root))
    return lines


def global_source_lines(run_root, n_org: int, scored: bool = True,
                        operon_db=None) -> list[str]:
    """Returns the footer source list for a run-level (pangenome) figure."""
    lines = []
    if scored:
        fname = _CONF_FINAL.split("/")[-1]
        lines.append(f"*/scoring/{fname}  (all {n_org} genomes)  [output]")
    if operon_db is not None:
        lines.append(rel_or_host(operon_db, run_root))
    return lines


# Column glossaries drawn by draw_method_note; each "\n" starts a new paragraph.
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
    """Returns the width of s at fontsize as a fraction of the figure width, cached per word."""
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
    """Draws a fully justified note below the plot content and registers it for the crop."""
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
        """Returns the figure-fraction width of s."""
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
        # Greedy word-wrap to the available width.
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
    """Draws a small footer listing the files behind the figure, below the lowest content."""
    if not sources:
        return
    run_root = Path(run_root)
    # This label's ":" is the only colon separator; everything else uses "|".
    header = (f"Files used to build this figure:  paths relative to output/{run_root.name}/  "
              f"|  [host] = reference file on the compute host")
    fw, fh = fig.get_figwidth(), max(fig.get_figheight(), 3.0)
    approx = max(60, int(fw / (0.011 * fontsize)))
    body = "   |   ".join(sources)
    lines = [header] + (textwrap.wrap(body, width=approx) or [body])
    # Places the footer below the lowest content; savefig's crop extends to include it.
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
    """Returns a manifest description for a fig*_*.tsv companion file."""
    fig = name.split("_", 1)[0]
    rest = name.split("_", 1)[1].rsplit(".", 1)[0].replace("_", " ") if "_" in name else ""
    return f"exact data table plotted in {fig}: {rest}".rstrip(": ")


def write_sources_manifest(outdir: Path, run_root: Path, organisms: list[str],
                           operon_db: Path,
                           occ_reference: Path = DEFAULT_OCC_REFERENCE) -> None:
    """Writes figure-sources.tsv listing every input file and fig*_*.tsv with path, status and role."""
    from datetime import datetime
    outdir, run_root = Path(outdir), Path(run_root)
    rows = []

    def add(role: str, path, description: str) -> None:
        """Appends one manifest row with existence, size and modification time."""
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
        add("input: scored genes",
            _resolve_organism_file(base, _CONF_FINAL) or base / _CONF_FINAL,
            f"[{org}] final per-gene confidence table (C1-C4, preliminary, "
            "operon-adjusted, final score, tier) with operon_id and member "
            "counts -- the values shown in the plots")
        add("input: gene coordinates",
            _resolve_organism_file(base, _LABELED) or base / _LABELED,
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


# ---- small numeric helpers ----
def pca_svd(X: np.ndarray):
    """Returns (scores, loadings, explained_variance_ratio) of standardized numpy SVD PCA."""
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
    """Returns the operon-size band label for a member count."""
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
    """Draws value_col boxplots per ordered bin with the bin means overlaid as a line.

    Returns (labels, ns, medians, means) for the companion TSV."""
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
    # Each bin gets its own colour.
    for i, patch in enumerate(bp["boxes"]):
        patch.set_facecolor(CATEGORICAL[i % len(CATEGORICAL)])
        patch.set_alpha(0.85); patch.set_edgecolor("white"); patch.set_linewidth(1.5)
    medians = [float(np.median(x)) for x in data]
    means = [float(np.mean(x)) for x in data]
    # Mean overlay: black line with white diamonds.
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


# ---- gene-track drawing and descriptor index ----
# Genes are drawn as numbered arrows; full descriptors are listed in a table or index.
def tag_letter(i: int) -> str:
    """Returns the 1-based number used to tag gene i in diagrams and tables."""
    return str(i + 1)


# Physical gene-glyph size in inches, converted to data units so every arrow looks the same.
_ARROW_W_IN = 0.74     # total glyph length (horizontal)
_ARROW_HLEN_IN = 0.16  # head length
_ARROW_BODY_IN = 0.18  # body thickness (vertical)
_ARROW_HEAD_IN = 0.36  # head thickness (vertical, wider than body)
_TRACK_YRANGE = 1.9    # data height of one arrow row (== row_pitch)
_TRACK_ROW_IN = 1.28   # inches per arrow row; the atlas driver uses the same value
# Inches per data-y-unit that pin_track_scale() pins every arrow track to after layout.
_TRACK_IN_PER_Y = _TRACK_ROW_IN / _TRACK_YRANGE


def _arrow_geometry(ax, span: int, nrows: int = 1):
    """Returns glyph geometry (hw, body_ht, head_ht, head_len) in data units.

    Height uses the pinned _TRACK_IN_PER_Y scale; width uses the measured axis width."""
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
    """Resets a track's ylim after layout so its inches-per-data-unit equals _TRACK_IN_PER_Y.

    Called once per track axis after the figure's final layout."""
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
    # Never frames tighter than the drawn rows.
    rng = max(rng, info["nrows"] * _TRACK_YRANGE)
    ax.set_ylim(info["y_top"] - rng, info["y_top"])


def member_descriptor(m) -> str:
    """Returns the cleaned descriptor prefixed with its annotation source, for drawing and sizing tables."""
    full = clean_descriptor(m.get("label", "")) or (m.get("label", "") or "—")
    src = (m.get("source") or "").strip()
    if src and full != "—":
        full = f"{src}: {full}"
    return full


def member_table_units(members, desc_wrap: int = None) -> float:
    """Returns the line units a gene table needs: header plus each wrapped source-prefixed descriptor."""
    w = desc_wrap or _TABLE_DESC_WRAP
    return _TABLE_HEADER_UNITS + sum(len(wrap_desc(member_descriptor(m), w)) for m in members)


def _gene_arrow(cx, y, strand, col, geom):
    """Returns a block-arrow Polygon pointing 5'->3' (right for +, left for -)."""
    hw, bt, ht, hl = geom
    if str(strand) == "-":
        tip, sh = cx - hw, cx - hw + hl
        pts = [(cx + hw, y + bt), (sh, y + bt), (sh, y + ht), (tip, y),
               (sh, y - ht), (sh, y - bt), (cx + hw, y - bt)]
    else:
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
    """Draws genes as numbered block arrows on a backbone with operon bands and gap labels.

    min_span fixes the slot count so arrows keep their size; per_row wraps long
    operons onto stacked rows; strand_rows is ignored. Returns (tag, descriptor,
    score, location, colour) tuples for render_gene_table()."""
    n = len(members)
    if n == 0:
        ax.axis("off")
        return []
    gene_col = [CATEGORICAL[(letter_offset + i) % len(CATEGORICAL)] for i in range(n)]
    # Wraps onto rows of at most per_row genes only when the operon is longer.
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
    # Fixed row height keeps the data range at nrows * row_pitch, so arrow size is constant.
    row_pitch = _TRACK_YRANGE

    def rc(k):
        """Returns (row, local column, baseline y) of gene k."""
        r = k // per
        return r, k - r * per, -r * row_pitch

    # Operon bands over same-operon runs within a row, drawn first so arrows sit on top.
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

    # Genome backbone behind the arrows, one segment per row.
    for r in range(nrows):
        y_r = -r * row_pitch
        if nrows > 1:
            last_col = (min((r + 1) * per, n) - 1) - r * per
            # Extends 0.60 on both sides of the row's genes.
            bx0, bx1 = off - 0.60, last_col + off + 0.60
        else:
            bx0, bx1 = x0, x1
        ax.plot([bx0, bx1], [y_r, y_r], color="#b7bdc6", lw=1.6, zorder=1,
                solid_capstyle="round")

    # Intergenic gaps (bp) under the backbone, within a row only.
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
    # Data range is exactly nrows * row_pitch to match _arrow_geometry's scale.
    y_top = ht + 0.48 + top_extra
    ax.set_ylim(y_top - nrows * row_pitch, y_top)
    # Stores what pin_track_scale() needs after layout.
    ax._track_scale = {"y_top": y_top, "nrows": nrows}
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ("top", "right", "bottom", "left"):
        ax.spines[s].set_visible(False)
    return entries


def _fmt_num(v):
    """Returns v as a float, or None for None/NaN."""
    return None if v is None or (isinstance(v, float) and np.isnan(v)) else float(v)


def _score_breakdown(score) -> str:
    """Returns a score summary such as 'C3 0.65   confidence = (prelim 0.78, operon +0.05, final 0.83)'.

    Missing parts are omitted; a plain float gives '(0.83)'."""
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
    """Draws the gene-tag index of full descriptors with score breakdowns.

    entries are (tag, descriptor, score); two_line puts the breakdown on a second line."""
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


_TABLE_DESC_WRAP = 46   # descriptor wrap width (chars)
_TABLE_HEADER_UNITS = 2.4   # two header lines (name + factor sub-line) plus padding


def wrap_desc(desc: str, width: int = _TABLE_DESC_WRAP) -> list[str]:
    """Returns the descriptor wrapped to width, or an em-dash line when empty."""
    return textwrap.wrap(desc or "—", width=width) or ["—"]


def table_line_units(entries, desc_wrap: int = _TABLE_DESC_WRAP) -> float:
    """Returns the line units a gene table needs, for sizing its axis before drawing."""
    return _TABLE_HEADER_UNITS + sum(len(wrap_desc(e[1], desc_wrap)) for e in entries)


def _conflict_code(score) -> str:
    """Returns conflict tags joined by '+': EC (C4 < 1), ambig (ambiguous operon), desc; else an em-dash."""
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
    """Returns a one-line review reason (at most 48 chars), or an em-dash when not flagged."""
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
    """Returns column x-positions (axis fractions) for the full-breakdown table between left and right.

    Fixed columns take their content width; descriptor and review reason share the
    rest. Computed once per page so every table on it aligns."""
    fs, fsn, fsf = fontsize, fontsize - 0.4, fontsize - 2.1
    Win = page_width_in or 18.0
    def _wf(s, pt):
        return (len(str(s)) * 0.60 * pt / 72.0) / Win
    GAP = 0.24 / Win
    def _numw(h, sub):
        return max(_wf(h, fs), _wf(sub, fsf), _wf("+0.00", fsn))
    def _numw2(h, sub):        # column holding an "adj/hyb" pair such as 0.00/0.00
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
    """Returns the largest wrap width (<= desc_wrap, >= 20) at which the descriptor column fits.

    Uses the same width model as breakdown_col_layout."""
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
    # Width left for the descriptor after the fixed and reason columns.
    w_desc = (right - left) - fixed - reason_nat
    char_w = 0.60 * fs / 72.0 / Win
    n = int(w_desc / char_w) if char_w > 0 else desc_wrap
    return max(20, min(desc_wrap, n))


def render_gene_table(ax, entries, fontsize: float = 7.6, show_location: bool = True,
                      show_scores: bool = True, desc_wrap: int = _TABLE_DESC_WRAP,
                      full_breakdown: bool = False, page_width_in: float = None,
                      col_layout: dict = None) -> None:
    """Draws a gene table: tag, wrapped descriptor, location, score columns and review flag.

    entries are draw_gene_track tuples; full_breakdown adds C1, C4, type and review reason."""
    ax.axis("off")
    if not entries:
        return
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    have_scores = show_scores and any(
        isinstance(e[2], dict) and _fmt_num(e[2].get("final")) is not None for e in entries)
    have_loc = show_location and any(len(e) > 3 and e[3] for e in entries)
    have_review = show_scores and any(
        isinstance(e[2], dict) and _fmt_num(e[2].get("prelim")) is not None
        and _fmt_num(e[2].get("operon")) is not None for e in entries)
    fs, fsn, fsf = fontsize, fontsize - 0.4, fontsize - 2.1
    wrapped = [wrap_desc(e[1], desc_wrap) for e in entries]

    # Column x-positions.
    if full_breakdown:
        # A page-level col_layout keeps every table on the page aligned.
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
        """Draws one table cell at x on the current row."""
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
    # Sub-header line naming what each column is a function of.
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
                # Non-CDS features (RNA, prophage) are shown in a warm tone.
                tcol = "#555555" if fty.lower() == "cds" else "#a86400"
                ax.text(x_type, y0, fty, ha="left", va="top", fontsize=fsn, color=tcol)
        if have_scores and isinstance(score, dict):
            c2 = _fmt_num(score.get("c2")); c3 = _fmt_num(score.get("c3"))
            p = _fmt_num(score.get("prelim"))
            op = _fmt_num(score.get("operon")); f = _fmt_num(score.get("final"))
            c1 = _fmt_num(score.get("c1")); c4 = _fmt_num(score.get("c4"))
            # Drivers below the 0.5 neutral line are coloured LOWERED.
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
                # C4 < 1 means an EC conflict was applied.
                ax.text(x_c4, y0, f"{c4:.2f}", va="top", fontsize=fsn,
                        color=(LOWERED if c4 < 0.999 else "#333333"))
            if p is not None:
                ax.text(x_pre, y0, f"{p:.2f}", va="top", fontsize=fsn, color="#333333")
            if p is not None and op is not None:
                if full_breakdown:
                    # Shows the raw boost C2*C3 (adj/hyb), not final - prelim, which hides the clip at 1.
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
                    d = op - p                        # effective delta in the compact table
                    ax.text(x_op, y0, f"{d:+.2f}", va="top", fontsize=fsn,
                            color=(RAISED if d > 0 else (LOWERED if d < 0 else "#333333")))
            if f is not None:
                fh = _fmt_num(score.get("final_hybrid"))
                txt = f"{f:.2f}/{fh:.2f}" if (full_breakdown and fh is not None) else f"{f:.2f}"
                ax.text(x_fin, y0, txt, va="top", fontsize=fsn,
                        fontweight="bold", color="#111111")
        if have_review and isinstance(score, dict):
            if full_breakdown:
                # Uses the scored needs_review flag plus a readable reason.
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
                    if (op_r - p_r) <= -0.1:  # operon context lowered the score by >= 0.1
                        ax.text(x_rev, y0, "yes", va="top", fontsize=fsn, fontweight="bold",
                                color=LOWERED)
                    else:
                        ax.text(x_rev, y0, "no", va="top", fontsize=fsn, color="#5a9e6f")
        y -= len(wlines) * line_h


# ---- operon pages, ported from the genome viewer's operon map ----
# Mirrors operonFigureSVG in viz/gen_genome_viewer.py, in viewer pixels (1 px = 1/96 in).
VIEWER_TIER_COL = ["#1F77FF", "#00B84D", "#FFCC00", "#FF8C00", "#EE2233"]
VIEWER_OPERON_CYCLE = ["#1F77FF", "#FF8C00", "#00B84D", "#B65CFF", "#00C2D1", "#EE2233"]
VIEWER_NONCODE = "#d5d5d5"
VIEWER_REVIEW = "#c0143c"
_VW = 1480                      # page width (px)
_VMX = 32                       # side margin (px)
_V_ROW = 21                     # table row height (px)
_V_SERIF = ["Times New Roman", "Times", "Nimbus Roman", "Liberation Serif", "DejaVu Serif"]


def _v_tier_index(tier) -> int:
    """Returns the tier's index in CONF_TIER_ORDER, or -1 (non-coding/unknown)."""
    t = (str(tier or "")).strip().lower()
    return CONF_TIER_ORDER.index(t) if t in CONF_TIER_ORDER else -1


def _v_review_short(reason) -> str:
    """Returns the short trigger tag(s) of a review reason, as the viewer's revShort does."""
    if not reason or str(reason).strip() in ("", "nan", "-"):
        return "yes"
    s, t = str(reason).lower(), []
    if "ec conflict" in s:
        t.append("EC conflict")
    if "low confidence" in s:
        t.append("low conf.")
    if "operon inference ambig" in s:
        t.append("operon ambig.")
    return "; ".join(t) if t else (str(reason) if len(str(reason)) <= 20 else str(reason)[:18] + "…")


def _v_num(v):
    """Returns v as a float, or None for missing/NaN/unparseable values."""
    try:
        f = float(v)
        return None if f != f else f
    except (TypeError, ValueError):
        return None


def _v_dec(v) -> str:
    """Returns v with two decimals, or an em-dash."""
    f = _v_num(v)
    return "—" if f is None else f"{f:.2f}"


def _v_flagged(m) -> bool:
    """Returns True when the member is flagged for review."""
    return str(m.get("needs_review") or "").strip().lower() in ("yes", "true", "1")


def operon_block_px(n_genes: int) -> int:
    """Returns the pixel height of one operon block (headings, arrows and table)."""
    return 104 + 30 + 48 + 24 + _V_ROW * max(n_genes, 1) + 40


def _v_block(ax, y0, members, heading, detail):
    """Draws one operon block at vertical offset y0 (px) and returns its height."""
    T = lambda x, y, text, px, bold=False, ha="left": ax.text(
        x, y, text, fontsize=px * 0.75, fontweight="bold" if bold else "normal",
        ha=ha, va="baseline", color="#000000", family=_V_SERIF)
    genes = [m for m in members if _v_num(m.get("start")) is not None and _v_num(m.get("end")) is not None]
    genes.sort(key=lambda m: min(_v_num(m["start"]), _v_num(m["end"])))
    n = len(genes)
    if not n:
        return 0
    s_ = lambda m: min(_v_num(m["start"]), _v_num(m["end"]))
    e_ = lambda m: max(_v_num(m["start"]), _v_num(m["end"]))
    lo, hi = min(s_(g) for g in genes), max(e_(g) for g in genes)
    rng = max(1.0, hi - lo)
    aw = _VW - 2 * _VMX
    sx = lambda v: _VMX + aw * (v - lo) / rng
    arr_y, arr_h = y0 + 104, 30
    head_y = arr_y + arr_h + 48

    T(_VMX, y0 + 30, heading, 20, bold=True)
    if detail:
        T(_VMX, y0 + 52, detail, 13)
    T(_VMX, y0 + 72, f"{n} genes  |  {(hi - lo) / 1000:.1f} kb region  |  arrow length ∝ gene length"
                     "  |  intergenic distances shown below arrows", 12.5)
    ax.plot([_VMX, _VW - _VMX], [arr_y + arr_h / 2] * 2, color="#cccccc", lw=1.2 * 0.75, zorder=1)

    spans, prev_end = [], None
    for k, g in enumerate(genes):
        x0, x1 = sx(s_(g)), sx(e_(g))
        if x1 - x0 < 7:
            x1 = x0 + 7
        if prev_end is not None and x0 < prev_end + 2:
            w = x1 - x0
            x0, x1 = prev_end + 2, prev_end + 2 + w
        prev_end = x1
        spans.append((x0, x1))
        hd = min(12.0, (x1 - x0) * 0.5)
        ti = _v_tier_index(g.get("tier"))
        col = VIEWER_NONCODE if ti < 0 else VIEWER_OPERON_CYCLE[k % len(VIEWER_OPERON_CYCLE)]
        yt, yb, ym = arr_y, arr_y + arr_h, arr_y + arr_h / 2
        plus = str(g.get("strand") or "+").strip() != "-"
        pts = ([(x0, yt), (x1 - hd, yt), (x1, ym), (x1 - hd, yb), (x0, yb)] if plus
               else [(x1, yt), (x0 + hd, yt), (x0, ym), (x0 + hd, yb), (x1, yb)])
        flag = _v_flagged(g)
        ax.add_patch(Polygon(pts, closed=True, facecolor=col,
                             edgecolor=VIEWER_REVIEW if flag else "#000000",
                             linewidth=(1.8 if flag else 0.7) * 0.75, joinstyle="miter", zorder=3))
        T((x0 + x1) / 2, arr_y - 10, str(k + 1), 11, ha="center")
    for k in range(n - 1):
        gap = int(s_(genes[k + 1]) - e_(genes[k]) - 1)
        lab = "‹1 bp" if gap <= 0 else (f"{gap} bp" if gap < 1000 else f"{gap / 1000:.1f} kb")
        T((spans[k][1] + spans[k + 1][0]) / 2, arr_y + arr_h + 15, lab, 10, ha="center")

    cx = {"loc": 432, "c1": 710, "c2": 770, "c3": 903, "c4": 965, "fin": 1090, "tier": 1180, "rev": _VW - _VMX}
    ax.plot([_VMX, _VW - _VMX], [head_y + 7] * 2, color="#000000", lw=0.8 * 0.75)
    T(_VMX, head_y, "#", 11.5, bold=True)
    T(_VMX + 21.5, head_y, "map", 8.5, bold=True, ha="center")
    T(_VMX + 41.5, head_y, "tier", 8.5, bold=True, ha="center")
    T(_VMX + 58, head_y, "gene product", 11.5, bold=True)
    T(cx["loc"], head_y, "location (bp)", 11.5, bold=True)
    for key, lab in (("c1", "C1"), ("c2", "C2"), ("c3", "C3 adj/hyb"), ("c4", "C4"),
                     ("fin", "final adj/hyb"), ("tier", "tier"), ("rev", "review")):
        T(cx[key], head_y, lab, 11.5, bold=True, ha="right")
    for k, g in enumerate(genes):
        y = head_y + 24 + _V_ROW * k
        ti = _v_tier_index(g.get("tier"))
        c = VIEWER_NONCODE if ti < 0 else VIEWER_OPERON_CYCLE[k % len(VIEWER_OPERON_CYCLE)]
        ct = VIEWER_NONCODE if ti < 0 else VIEWER_TIER_COL[ti]
        tn = "non-coding" if ti < 0 else CONF_TIER_ORDER[ti]
        T(_VMX, y, str(k + 1), 11)
        for dx, fill in ((16, c), (36, ct)):
            ax.add_patch(Rectangle((_VMX + dx, y - 9), 11, 11, facecolor=fill,
                                   edgecolor="#000000", linewidth=0.5 * 0.75, zorder=3))
        name = str(g.get("label") or "").strip() or "(unnamed)"
        T(_VMX + 58, y, name[:54], 12)
        plus = str(g.get("strand") or "+").strip() != "-"
        T(cx["loc"], y, f"{int(s_(g)):,}–{int(e_(g)):,} {'+' if plus else '−'}", 11)
        T(cx["c1"], y, _v_dec(g.get("c1")), 12, ha="right")
        T(cx["c2"], y, _v_dec(g.get("c2")), 12, ha="right")
        T(cx["c3"], y, f"{_v_dec(g.get('c3'))}/{_v_dec(g.get('c3_hybrid'))}", 12, ha="right")
        T(cx["c4"], y, _v_dec(g.get("c4")), 12, ha="right")
        T(cx["fin"], y, f"{_v_dec(g.get('operon_adjusted'))}/{_v_dec(g.get('operon_adjusted_hybrid'))}",
          12, bold=True, ha="right")
        T(cx["tier"], y, tn, 11, ha="right")
        T(cx["rev"], y, _v_review_short(g.get("review_reason")) if _v_flagged(g) else "", 11, ha="right")
    return operon_block_px(n)


def render_operon_page(outpath, blocks, *, org_label, suptitle, dpi=288, **_legacy) -> None:
    """Writes one page of operon blocks under a title and organism, with the tier key at the foot.

    blocks are {"members", "heading" or "title", "detail"} dicts; extra keyword arguments are ignored."""
    blocks = [b for b in blocks if b.get("members")]
    if not blocks:
        return
    top = 30 + (22 if org_label else 0) + 26
    body = sum(operon_block_px(len(b["members"])) for b in blocks)
    H = top + body + 40
    fig = plt.figure(figsize=(_VW / 96.0, H / 96.0))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, _VW)
    ax.set_ylim(H, 0)
    ax.axis("off")
    ax.text(_VMX, 30, suptitle, fontsize=22 * 0.75, fontweight="bold", va="baseline", family=_V_SERIF)
    if org_label:
        ax.text(_VMX, 52, str(org_label).replace("_", " "), fontsize=14 * 0.75, va="baseline",
                style="italic", fontweight="normal", family=_V_SERIF)
    y = top
    for b in blocks:
        heading = b.get("heading") or b.get("title") or (b["members"][0].get("operon_id") or "")
        y += _v_block(ax, y, b["members"], heading, b.get("detail", ""))
    lx, ly = _VMX, H - 26
    for i, nm in enumerate(CONF_TIER_ORDER + ["non-coding"]):
        ax.add_patch(Rectangle((lx, ly), 12, 12, facecolor=VIEWER_TIER_COL[i] if i < 5 else VIEWER_NONCODE,
                               edgecolor="#000000", linewidth=0.5 * 0.75))
        ax.text(lx + 16, ly + 10, nm, fontsize=11 * 0.75, va="baseline", fontweight="normal", family=_V_SERIF)
        lx += 16 + len(nm) * 6.3 + 20
    ax.text(lx + 4, ly + 10, "red outline = flagged for review", fontsize=11 * 0.75, va="baseline",
            fontweight="normal", family=_V_SERIF)
    path = Path(outpath)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, facecolor="white")
    plt.close(fig)
    print(f"[reportfig] wrote {path.name}", file=sys.stderr)


def operon_members_informative(members_in_order: str, min_informative: int | None = None):
    """Returns True for an operon of >=2 members that are all (or >= min_informative) informative.

    Uses scoring's is_uninformative gate; "-" members count as uninformative."""
    members = [m.strip() for m in (members_in_order or "").split(" -> ") if m.strip()]
    if len(members) < 2:
        return False
    informative = [(m not in ("-", "")) and not is_uninformative(m) for m in members]
    if min_informative is None:
        return all(informative)
    return sum(informative) >= min_informative


def operon_to_members(op_row) -> list[dict]:
    """Returns draw_gene_track member dicts, with per-gene scores, from a build_operons() row."""
    out = []
    labels = op_row["member_labels"]

    def col(name):
        """Returns the named list column, or Nones when absent."""
        return op_row[name] if name in op_row else [None] * len(labels)

    confs, prelims = col("confidences"), col("preliminaries")
    opadj, c3s, srcs = col("operon_adjusteds"), col("c3s"), col("sources")
    opadj_h, c3s_h = col("operon_adjusteds_hybrid"), col("c3s_hybrid")
    c2s = col("c2s")
    c1s, c4s = col("c1s"), col("c4s")
    nrs = col("needs_reviews")
    ftys = col("feature_types")
    reasons = col("review_reasons")
    tiers = col("tiers")
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
            "tier": tiers[i] if i < len(tiers) else None,
        })
    return out
