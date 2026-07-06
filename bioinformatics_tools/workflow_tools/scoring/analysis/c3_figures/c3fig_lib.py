#!/usr/bin/env python3
"""c3fig_lib.py — shared helpers for the C3 comprehensive figure suite.

Provides:
  * is_uninformative()  — EXACT copy of the authoritative gate in
    scoring/score-c1-tool-coverage.py + labeling/assign-canonical-label.py.
    A gene is "uninformative" (hypothetical / unknown-function) iff the
    functional text of its best_consensus_product_descriptor matches this
    gate; otherwise it is "informative" (carries a real function name).
  * clean_descriptor()  — strip the leading "SOURCE: " tag and the
    "raw ## human" duplication so genes group by functional identity.
  * data loading         — JOIN labeled-genes.tsv (coordinates, descriptor,
    aa_seq) with labeled-genes-operon-info.tsv (operon_id, member_count,
    position, probability) on feature_id, across all organisms.
  * matplotlib styling    — Times-New-Roman-compatible serif, everything bold.
  * pca()                — numpy-SVD PCA (sklearn is unavailable on py3.6).

The cache built by c3fig_00_build_cache.py is a pickle of a pandas DataFrame
(one row per protein-coding gene, all organisms) so every figure script loads
the parsed data in <1 s instead of re-reading ~2 GB of TSV.
"""
import csv
import hashlib
import pickle
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

csv.field_size_limit(10_000_000)

# ─────────────────────────────────────────────────────────────────────────────
# UNINFORMATIVE HIT CATEGORY — kept byte-for-byte identical to the gate in
# scoring/score-c1-tool-coverage.py (which mirrors labeling/assign-canonical-
# label.py). If you edit one, edit all three.
# ─────────────────────────────────────────────────────────────────────────────
_UNINFORMATIVE = frozenset({
    "", "-", ".", "na", "n/a", "none", "null",
    "unknown", "uncharacterized", "uncharacterised", "putative", "predicted",
    "hypothetical protein", "conserved hypothetical protein", "conserved protein",
    "conserved domain protein", "predicted protein", "function unknown",
    "domain of unknown function", "general function prediction only",
    "poorly characterized", "open reading frame",
})
_UNINFORMATIVE_PREFIXES = (
    "domain of unknown function", "protein of unknown function",
    "family of unknown function", "region of unknown function",
    "repeat of unknown function", "module of unknown function",
    "duf", "upf", "uncharacteri", "putative uncharacteri",
    "conserved hypothetical", "hypothetical", "unknown protein",
    "unknown function", "orf", "pfam uncharacteri",
)
_UNKNOWN_FUNC_RE = re.compile(
    r'^\s*(?:(?:bacterial|viral|archaeal|eukaryotic|fungal|plant|marine|'
    r'transmembrane|integral membrane|membrane)\s+)?'
    r'(?:domain|protein|family|repeat|region|module)\s+of\s+unknown\s+function',
    re.IGNORECASE,
)
_DB_ID_HYPOTHETICAL_RE = re.compile(
    r'^(?:fig\d+|tigr\d+)[:\s].*(?:hypothetical|conserved hypothetical)', re.IGNORECASE,
)
_UPF_ONLY_RE = re.compile(r'^\s*belongs to the upf\d+', re.IGNORECASE)
_BARE_DUF_RE = re.compile(r'^\s*(?:pfam:)?\(?duf\d+\)?(?:\s+(?:family|domain))?\s*$', re.IGNORECASE)
_PROTEIN_DOMAINS_DUF_RE = re.compile(r'^\s*protein containing domains?\s+duf', re.IGNORECASE)
_HYPOTHETICAL_ANYWHERE_RE = re.compile(r'\bhypothetical\b', re.IGNORECASE)
_PROTEIN_CONSERVED_IN_BACTERIA_RE = re.compile(r'\bprotein\s+conserved\s+in\s+bacteria\b', re.IGNORECASE)
_INTEGRAL_MEMBRANE_PROTEIN_RE = re.compile(r'^\s*integral\s+membrane\s+protein\s*$', re.IGNORECASE)
_FUNCTION_SIGNAL_RE = re.compile(
    r'(ase\b|transport|permease|pump|export|import|channel|carrier|symport|antiport|'
    r'bind|synth|kinas|reductas|hydrolas|transferas|isomeras|ligas|lyas|mutas|oxidas|'
    r'dehydrogen|regulat|repressor|activator|\bfactor\b|receptor|sensor|subunit|ribosom|'
    r'polymeras|oxidoreduct|enzyme|proteas|peptidas|nucleas|phosphatas|efflux|resistance|'
    r'virulence|toxin|flippase|homeostasis|assembly|motility|adhesin|chaperone|'
    r'helicas|topoisomeras|gyrase|recombinas|integras|transposas|methylas|glycosyl|'
    r'dismutas|catalas|peroxidas|cytochrome|ferredoxin|cytoskelet|cell division|'
    r'degradation|tolerance|translation|utilization)', re.IGNORECASE,
)
_UNCHARACTERIZED_RE = re.compile(r'\buncharacteri[sz]ed\b', re.IGNORECASE)
_LOCUS_OR_TAG_RE = re.compile(r'^\(?(?:duf\d+|upf\d+|[a-z]{1,5}\d{0,4}[a-z]?\d{0,4})\)?$', re.IGNORECASE)


def is_uninformative(val: str) -> bool:
    """True iff ``val`` falls in the UNINFORMATIVE HIT CATEGORY."""
    v = (val or "").strip().lower()
    if not v or v in _UNINFORMATIVE:
        return True
    for prefix in _UNINFORMATIVE_PREFIXES:
        if v.startswith(prefix):
            return True
    if (_UNKNOWN_FUNC_RE.match(v) or _DB_ID_HYPOTHETICAL_RE.match(v)
            or _UPF_ONLY_RE.match(v) or _BARE_DUF_RE.match(v)
            or _HYPOTHETICAL_ANYWHERE_RE.search(v)
            or _PROTEIN_CONSERVED_IN_BACTERIA_RE.search(v)
            or _INTEGRAL_MEMBRANE_PROTEIN_RE.match(v)
            or _PROTEIN_DOMAINS_DUF_RE.match(v)):
        return True
    if not _FUNCTION_SIGNAL_RE.search(v):
        if "of unknown function" in v:
            return True
        if _UNCHARACTERIZED_RE.search(v):
            return True
        if v.startswith("conserved protein"):
            rest = v[len("conserved protein"):].strip(" ,.;:-")
            if not rest or _BARE_DUF_RE.match(rest) or _LOCUS_OR_TAG_RE.match(rest):
                return True
    return False


# ─────────────────────────────────────────────────────────────────────────────
# Descriptor cleaning
# ─────────────────────────────────────────────────────────────────────────────
# Leading "SOURCE: " tags that the labeling pipeline prepends to the functional
# text.  We strip these so the SAME function annotated via different evidence
# sources groups together, and so is_uninformative() sees the functional text.
_SOURCE_PREFIX_RE = re.compile(
    r'^(?:'
    r'JCVI|NCBIFAM|NCBI Protein Cluster \(PRK\)|PRK|TIGR|TIGRFAM|PGAP|HAMAP|'
    r'PIRSF|UniProt|UNIPROT|KEGG|EGGNOG|eggNOG|PFAM|Pfam|CDD|COG|RAST|MEROPS|'
    r'TCDB|DBCAN|SwissProt|FIG\d+'
    r')\s*:\s*',
    re.IGNORECASE,
)


def clean_descriptor(desc: str) -> str:
    """Return the functional text: strip a leading SOURCE tag and collapse the
    "raw ## human-readable" duplication (keep the human-readable side)."""
    d = (desc or "").strip()
    if not d:
        return ""
    # "5S rRNA ## 5S ribosomal RNA" -> keep the more descriptive right side
    if " ## " in d:
        left, right = d.split(" ## ", 1)
        d = right.strip() if right.strip() else left.strip()
    # strip a single leading SOURCE: tag (may repeat, e.g. "JCVI: PRK: x")
    prev = None
    while prev != d:
        prev = d
        d = _SOURCE_PREFIX_RE.sub("", d).strip()
    return d


def sha256_hash(seq: str) -> str:
    return hashlib.sha256((seq or "").strip().upper().encode("utf-8")).hexdigest()


# ─────────────────────────────────────────────────────────────────────────────
# Data loading
# ─────────────────────────────────────────────────────────────────────────────
NON_OPERON_TOKENS = {"NOT_IN_AN_OPERON", "NOT_APPLICABLE_NON_CODING", "UNKNOWN", ""}


def discover_organisms(labeling_root: Path):
    """Yield (organism, labeled_path, operon_path) for every organism that has
    both a labeled-genes.tsv and its operon-info sibling."""
    out = []
    for labeled_path in sorted(labeling_root.glob("**/labeling/labeled-genes.tsv")):
        operon_path = labeled_path.with_name("labeled-genes-operon-info.tsv")
        if operon_path.is_file():
            organism = labeled_path.parent.parent.name
            out.append((organism, labeled_path, operon_path))
    return out


def _to_int(v, default=0):
    try:
        return int(float(str(v).strip()))
    except (ValueError, TypeError):
        return default


def _to_float(v, default=float("nan")):
    try:
        return float(str(v).strip())
    except (ValueError, TypeError):
        return default


def load_organism(organism: str, labeled_path: Path, operon_path: Path,
                  compute_hash: bool = True) -> pd.DataFrame:
    """JOIN the two per-organism TSVs on feature_id -> one row per gene."""
    # Pass 1: operon-info (small, authoritative for operon membership)
    operon = {}
    with open(operon_path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            fid = (row.get("feature_id") or "").strip()
            if not fid:
                continue
            operon[fid] = {
                "operon_id": (row.get("operon_id") or "").strip(),
                "member_count": _to_int(row.get("operon_member_count")),
                "position_in_operon": _to_int(row.get("operon_gene_position_in_operon")),
                "operon_prob": _to_float(row.get("operon_probability_geometric_mean")),
            }

    # Pass 2: labeled-genes (coordinates, descriptor, aa_seq)
    records = []
    with open(labeled_path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            fid = (row.get("feature_id") or "").strip()
            if not fid:
                continue
            full_desc = (row.get("best_consensus_product_descriptor") or "").strip()
            clean = clean_descriptor(full_desc)
            aa_seq = (row.get("aa_seq") or "").strip()
            op = operon.get(fid, {})
            operon_id = op.get("operon_id", "")
            in_operon = operon_id.startswith("operon_")
            rec = {
                "organism": organism,
                "feature_id": fid,
                "full_descriptor": full_desc,
                "clean_descriptor": clean,
                "aa_hash": sha256_hash(aa_seq) if (compute_hash and aa_seq) else "",
                "aa_length": _to_int(row.get("aa_length")),
                "start": _to_int(row.get("gene_start")),
                "end": _to_int(row.get("gene_end")),
                "strand": (row.get("RAST_strand") or "").strip() or "+",
                "operon_id": operon_id if in_operon else "",
                "in_operon": in_operon,
                "member_count": op.get("member_count", 0),
                "position_in_operon": op.get("position_in_operon", 0),
                "operon_prob": op.get("operon_prob", float("nan")),
                "uninformative": is_uninformative(clean),
            }
            records.append(rec)
    return pd.DataFrame.from_records(records)


def build_gene_table(labeling_root: Path, compute_hash: bool = True) -> pd.DataFrame:
    frames = []
    for organism, lp, op in discover_organisms(labeling_root):
        print(f"[c3fig_lib] loading {organism} ...", file=sys.stderr)
        frames.append(load_organism(organism, lp, op, compute_hash=compute_hash))
    if not frames:
        raise SystemExit("[c3fig_lib] no organisms found under " + str(labeling_root))
    df = pd.concat(frames, ignore_index=True)
    df = attach_component_scores(df, labeling_root)
    return df


# ─────────────────────────────────────────────────────────────────────────────
# Deterministic confidence components joined per gene from
# scored-labeled-genes-confidence-final.tsv:
#   C1 = database/tool coverage (informative annotation sources / 7)
#   C2 = operon geometric-mean probability (0.5 = neutral for non-operonic genes)
#   C4 = EC-number agreement across tools
# C3 (operonic co-occurrence context) is DELIBERATELY EXCLUDED from this join:
# C3 is itself derived from the cross-organism co-occurrence signal that this
# whole suite analyses, so relating C3 to our operon co-occurrence results would
# be circular.  C1/C2/C4 are the deterministic components we compare against the
# operon structure derived elsewhere in the suite.
# ─────────────────────────────────────────────────────────────────────────────
COMPONENT_COLS = ["c1_score", "c2_score", "c4_score"]
COMPONENT_LABELS = {
    "c1_score": "C1 database coverage",
    "c2_score": "C2 operon probability",
    "c4_score": "C4 EC agreement",
}
OPERON_COMPONENT_COLS = ["mean_c1", "mean_c2", "mean_c4"]
_COMPONENT_SRC = {
    "c1_score": "c1_score",
    "c2_score": "c2_score_from_operon_probability",
    "c4_score": "c4_score",
}


def score_file_for(labeled_path: Path) -> Path:
    """<run>/<organism>/labeling/labeled-genes.tsv ->
       <run>/<organism>/scoring/scored-labeled-genes-confidence-final.tsv"""
    return (labeled_path.parent.parent / "scoring" /
            "scored-labeled-genes-confidence-final.tsv")


def load_component_scores(labeling_root: Path) -> pd.DataFrame:
    """Per-gene C1/C2/C4 for every organism (lean by-index column read)."""
    frames = []
    for organism, labeled_path, _ in discover_organisms(labeling_root):
        sp = score_file_for(labeled_path)
        if not sp.is_file():
            print(f"[c3fig_lib] no score file for {organism}", file=sys.stderr)
            continue
        recs = []
        with open(sp, newline="") as fh:
            rd = csv.reader(fh, delimiter="\t")
            hdr = next(rd)
            if "feature_id" not in hdr:
                continue
            i_fid = hdr.index("feature_id")
            idx = {k: hdr.index(v) for k, v in _COMPONENT_SRC.items() if v in hdr}
            need = max([i_fid] + list(idx.values()))
            for r in rd:
                if len(r) <= need:
                    continue
                rec = {"organism": organism, "feature_id": r[i_fid].strip()}
                for k, i in idx.items():
                    rec[k] = _to_float(r[i])
                recs.append(rec)
        frames.append(pd.DataFrame.from_records(recs))
    if not frames:
        return pd.DataFrame(columns=["organism", "feature_id"] + COMPONENT_COLS)
    return pd.concat(frames, ignore_index=True)


def attach_component_scores(genes: pd.DataFrame, labeling_root: Path) -> pd.DataFrame:
    """LEFT-join C1/C2/C4 onto the gene table by (organism, feature_id)."""
    sc = load_component_scores(labeling_root)
    if sc.empty:
        for c in COMPONENT_COLS:
            genes[c] = float("nan")
        return genes
    return genes.merge(sc, on=["organism", "feature_id"], how="left")


def save_cache(df: pd.DataFrame, cache_path: Path):
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    with open(cache_path, "wb") as fh:
        pickle.dump(df, fh, protocol=pickle.HIGHEST_PROTOCOL)


def load_cache(cache_path: Path) -> pd.DataFrame:
    with open(cache_path, "rb") as fh:
        return pickle.load(fh)


# ─────────────────────────────────────────────────────────────────────────────
# Operon-level table (derived from the gene table)
# ─────────────────────────────────────────────────────────────────────────────
def build_operon_table(genes: pd.DataFrame) -> pd.DataFrame:
    """One row per operon (organism + operon_id) with composition/geometry."""
    op = genes[genes["in_operon"]].copy()
    rows = []
    for (organism, operon_id), g in op.groupby(["organism", "operon_id"]):
        g = g.sort_values("start")
        n = len(g)
        n_unf = int(g["uninformative"].sum())
        n_inf = n - n_unf
        starts = g["start"].to_numpy()
        ends = g["end"].to_numpy()
        strands = g["strand"].to_numpy()
        # intergenic gaps between consecutive genes (bp); negative -> overlap
        gaps = []
        for i in range(n - 1):
            gaps.append(int(starts[i + 1] - ends[i]))
        frac_unf = n_unf / n if n else 0.0
        if n_unf == 0:
            comp = "all_informative"
        elif n_inf == 0:
            comp = "all_hypothetical"
        elif n_inf > n_unf:
            comp = "majority_informative"
        elif n_unf > n_inf:
            comp = "majority_hypothetical"
        else:
            comp = "equal"
        n_plus = int((strands == "+").sum())
        n_minus = int((strands == "-").sum())
        same_strand = (n_plus == n) or (n_minus == n)
        rows.append({
            "organism": organism,
            "operon_id": operon_id,
            "size": n,
            "n_informative": n_inf,
            "n_uninformative": n_unf,
            "frac_uninformative": frac_unf,
            "composition": comp,
            "span_bp": int(ends.max() - starts.min()),
            "mean_gap_bp": float(np.mean(gaps)) if gaps else float("nan"),
            "median_gap_bp": float(np.median(gaps)) if gaps else float("nan"),
            "min_gap_bp": int(np.min(gaps)) if gaps else 0,
            "max_gap_bp": int(np.max(gaps)) if gaps else 0,
            "n_plus": n_plus,
            "n_minus": n_minus,
            "same_strand": same_strand,
            "dominant_strand": "+" if n_plus >= n_minus else "-",
            "mean_operon_prob": float(g["operon_prob"].mean()),
            "mean_aa_length": float(g["aa_length"].mean()),
        })
    op_table = pd.DataFrame.from_records(rows)
    # attach per-operon mean deterministic components (C1/C2/C4) if present
    have = [c for c in COMPONENT_COLS if c in op.columns]
    if have and len(op_table):
        means = (op.groupby(["organism", "operon_id"])[have].mean()
                   .rename(columns={c: "mean_" + c.replace("_score", "")
                                    for c in have}).reset_index())
        op_table = op_table.merge(means, on=["organism", "operon_id"], how="left")
    return op_table


def genome_sizes(genes: pd.DataFrame) -> pd.Series:
    """Approx genome size per organism = max gene end coordinate (bp)."""
    return genes.groupby("organism")["end"].max()


def build_adjacent_pairs(genes: pd.DataFrame) -> pd.DataFrame:
    """One row per consecutive (adjacent) gene pair inside an operon, with the
    intergenic gap, strand pattern and informativeness class. Used by the
    intergenic-distance, strand and adjacency figures."""
    op = genes[genes["in_operon"]]
    rows = []
    for (organism, operon_id), g in op.groupby(["organism", "operon_id"]):
        g = g.sort_values("start")
        starts = g["start"].to_numpy()
        ends = g["end"].to_numpy()
        strands = g["strand"].to_numpy()
        unf = g["uninformative"].to_numpy()
        n = len(g)
        for i in range(n - 1):
            gap = int(starts[i + 1] - ends[i])
            u1, u2 = bool(unf[i]), bool(unf[i + 1])
            if not u1 and not u2:
                cls = "info_info"
            elif u1 and u2:
                cls = "hypo_hypo"
            else:
                cls = "info_hypo"
            rows.append({
                "organism": organism,
                "operon_id": operon_id,
                "operon_size": n,
                "gap_bp": gap,
                "strand_left": strands[i],
                "strand_right": strands[i + 1],
                "strand_pattern": f"{strands[i]}/{strands[i + 1]}",
                "same_strand": strands[i] == strands[i + 1],
                "pair_class": cls,
            })
    return pd.DataFrame.from_records(rows)


# ─────────────────────────────────────────────────────────────────────────────
# Co-occurrence of genes (by product descriptor) inside operons.
#
# A "gene" here is identified by its clean_descriptor (functional name). Two
# genes CO-OCCUR when they belong to the same operon (order / distance / adjacency
# do NOT matter). Because hashes (exact aa sequences) essentially never recur
# across organisms, co-occurrence is measured on functional descriptors so it can
# be compared across the 21 organisms.
#
# Uninformative (hypothetical) descriptors are EXCLUDED from co-occurrence by
# default: "hypothetical protein" is not a gene identity — thousands of unrelated
# genes share it — so pairing on it is meaningless and would swamp real pairs.
# (Uninformative genes are still kept and counted everywhere else in the suite.)
# ─────────────────────────────────────────────────────────────────────────────
from collections import defaultdict as _defaultdict
from itertools import combinations as _combinations


def _operon_descriptor_sets(genes, informative_only=True):
    """Yield (organism, operon_id, sorted-unique-descriptors) for every operon."""
    op = genes[genes["in_operon"]]
    if informative_only:
        op = op[~op["uninformative"]]
    for (org, oid), g in op.groupby(["organism", "operon_id"]):
        descs = sorted({d for d in g["clean_descriptor"] if d})
        yield org, oid, descs


def build_cooccurrence_pairs(genes, informative_only=True):
    """Unordered co-occurring descriptor pairs within operons (any position).

    Returns DataFrame [desc_a, desc_b, n_operons, n_organisms, organisms] sorted
    by cross-organism spread then operon frequency. n_operons = number of operons
    (across all organisms) containing BOTH genes; n_organisms = number of distinct
    organisms in which the pair co-occurs in at least one operon."""
    op_count = _defaultdict(int)
    org_sets = _defaultdict(set)
    for org, _oid, descs in _operon_descriptor_sets(genes, informative_only):
        for a, b in _combinations(descs, 2):
            op_count[(a, b)] += 1
            org_sets[(a, b)].add(org)
    rows = []
    for (a, b), c in op_count.items():
        orgs = org_sets[(a, b)]
        rows.append({"desc_a": a, "desc_b": b, "n_operons": c,
                     "n_organisms": len(orgs),
                     "organisms": ";".join(sorted(short_label(o) for o in orgs))})
    df = pd.DataFrame.from_records(rows)
    if len(df):
        df = df.sort_values(["n_organisms", "n_operons"],
                            ascending=False).reset_index(drop=True)
    return df


def build_adjacent_cooccurrence(genes, informative_only=True):
    """Co-occurring descriptor pairs that are IMMEDIATELY ADJACENT in an operon.

    Returns DataFrame [desc_a, desc_b, n_adjacent, n_organisms, median_gap_bp,
    mean_gap_bp, frac_plus_strand]. Identical-descriptor neighbours (tandem
    duplicates) are skipped so each row is a pair of two distinct genes."""
    op = genes[genes["in_operon"]]
    stats = {}
    for (_org, _oid), g in op.groupby(["organism", "operon_id"]):
        g = g.sort_values("start")
        d = g["clean_descriptor"].to_numpy()
        s = g["start"].to_numpy(); e = g["end"].to_numpy()
        st = g["strand"].to_numpy(); un = g["uninformative"].to_numpy()
        for i in range(len(g) - 1):
            if informative_only and (un[i] or un[i + 1]):
                continue
            da, db = d[i], d[i + 1]
            if not da or not db or da == db:
                continue
            key = (da, db) if da < db else (db, da)
            rec = stats.setdefault(key, {"n": 0, "gaps": [], "orgs": set(),
                                         "plus": 0})
            rec["n"] += 1
            rec["gaps"].append(int(s[i + 1] - e[i]))
            rec["orgs"].add(_org)
            if st[i] == "+":
                rec["plus"] += 1
    rows = []
    for (a, b), r in stats.items():
        rows.append({"desc_a": a, "desc_b": b, "n_adjacent": r["n"],
                     "n_organisms": len(r["orgs"]),
                     "median_gap_bp": float(np.median(r["gaps"])),
                     "mean_gap_bp": float(np.mean(r["gaps"])),
                     "frac_plus_strand": r["plus"] / r["n"]})
    df = pd.DataFrame.from_records(rows)
    if len(df):
        df = df.sort_values(["n_organisms", "n_adjacent"],
                            ascending=False).reset_index(drop=True)
    return df


def build_kmember_sets(genes, kmax=4, informative_only=True):
    """Recurring co-occurring gene SETS of size k (k=2..kmax) within operons.

    For every operon we enumerate all k-subsets of its unique informative
    descriptors and count how many operons / organisms contain each subset.
    Returns DataFrame [k, members, n_operons, n_organisms] (members = ' + '
    joined sorted descriptors), sorted by k then cross-organism spread. Only
    subsets that recur (n_operons >= 2) are kept to bound the table size."""
    per_k_count = {k: _defaultdict(int) for k in range(2, kmax + 1)}
    per_k_orgs = {k: _defaultdict(set) for k in range(2, kmax + 1)}
    for org, _oid, descs in _operon_descriptor_sets(genes, informative_only):
        n = len(descs)
        for k in range(2, min(kmax, n) + 1):
            for combo in _combinations(descs, k):
                per_k_count[k][combo] += 1
                per_k_orgs[k][combo].add(org)
    rows = []
    for k in range(2, kmax + 1):
        for combo, c in per_k_count[k].items():
            if c < 2:
                continue
            rows.append({"k": k, "members": " + ".join(combo), "n_operons": c,
                         "n_organisms": len(per_k_orgs[k][combo])})
    df = pd.DataFrame.from_records(rows)
    if len(df):
        df = df.sort_values(["k", "n_organisms", "n_operons"],
                            ascending=[True, False, False]).reset_index(drop=True)
    return df


def load_or_build(name, builder, cache_dir):
    """Lazily load _cache/<name>.pkl, building + caching it on first use."""
    p = Path(cache_dir) / f"{name}.pkl"
    if p.is_file():
        return load_cache(p)
    df = builder()
    save_cache(df, p)
    return df


# ─────────────────────────────────────────────────────────────────────────────
# Operon feature matrix for PCA / clustering (Theme 7). One row per operon with
# the interpretable numeric features that could drive clustering: operon size,
# span, intergenic spacing, composition (fraction hypothetical), gene length,
# operon probability and host genome size. Strand homogeneity is dropped (every
# operon is single-strand -> zero variance). Skewed positive features are
# log-scaled; the signed intergenic gap uses a sign-preserving log.
# ─────────────────────────────────────────────────────────────────────────────
PCA_FEATURES = [
    ("size", "Operon size (genes)"),
    ("log_span_bp", "Operon span (log10 bp)"),
    ("signed_log_gap", "Mean intergenic gap (signed log)"),
    ("frac_uninformative", "Fraction hypothetical"),
    ("mean_aa_length", "Mean protein length (aa)"),
    ("mean_operon_prob", "Operon probability"),
    ("log_genome_size", "Genome size (log10 bp)"),
]


def operon_feature_matrix(operons, genes):
    """Return (F, feature_cols, feature_labels). F is a per-operon DataFrame with
    the transformed PCA features plus metadata columns (organism, size,
    frac_uninformative, composition) for colouring."""
    gsize = genome_sizes(genes)
    F = operons.copy()
    F["genome_size"] = F["organism"].map(gsize)
    F["log_span_bp"] = np.log10(F["span_bp"].clip(lower=1))
    F["log_genome_size"] = np.log10(F["genome_size"].clip(lower=1))
    F["signed_log_gap"] = np.sign(F["mean_gap_bp"]) * np.log1p(F["mean_gap_bp"].abs())
    cols = [c for c, _ in PCA_FEATURES]
    labels = [lab for _, lab in PCA_FEATURES]
    F = F.dropna(subset=cols).reset_index(drop=True)
    return F, cols, labels


# ─────────────────────────────────────────────────────────────────────────────
# Numpy-SVD PCA (sklearn unavailable on this py3.6)
# ─────────────────────────────────────────────────────────────────────────────
def pca(X: np.ndarray, n_components: int = 2):
    """Standardise columns, run SVD PCA.
    Returns (scores[n,k], loadings[features,k], explained_variance_ratio[k])."""
    X = np.asarray(X, dtype=float)
    mu = np.nanmean(X, axis=0)
    sd = np.nanstd(X, axis=0)
    sd[sd == 0] = 1.0
    Xs = (X - mu) / sd
    Xs = np.nan_to_num(Xs, nan=0.0)
    U, S, Vt = np.linalg.svd(Xs, full_matrices=False)
    scores = U * S
    ev = (S ** 2) / (Xs.shape[0] - 1)
    evr = ev / ev.sum()
    k = min(n_components, Vt.shape[0])
    return scores[:, :k], Vt[:k].T, evr[:k]


def response_surface(x, y, z, c=None, nbin=14, nx=44, ny=44,
                     min_count=10, smooth=1.0):
    """Build a smooth 3-D response surface z = f(x, y) from scattered data.

    Bins (x, y) into nbin*nbin cells and takes the MEAN z (and mean c) per cell.
    Cells with fewer than `min_count` samples are discarded so sparse cells cannot
    create spurious spikes; the surviving cell means are interpolated onto an
    nx*ny mesh and lightly Gaussian-smoothed (`smooth` = sigma in mesh cells).
    Optional 4th dimension `c` is returned as a matching mesh for colour mapping.
    Returns (Xi, Yi, Zi) or (Xi, Yi, Zi, Ci)."""
    from scipy.stats import binned_statistic_2d
    from scipy.interpolate import griddata
    from scipy.ndimage import gaussian_filter
    x = np.asarray(x, float); y = np.asarray(y, float); z = np.asarray(z, float)
    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    if c is not None:
        c = np.asarray(c, float); m = m & np.isfinite(c)
    x, y, z = x[m], y[m], z[m]
    if c is not None:
        c = c[m]

    def _fill(vals_valid, pts, Xi, Yi):
        Vi = griddata(pts, vals_valid, (Xi, Yi), method="linear")
        if np.isnan(Vi).any():
            Vn = griddata(pts, vals_valid, (Xi, Yi), method="nearest")
            Vi = np.where(np.isnan(Vi), Vn, Vi)
        if smooth:
            Vi = gaussian_filter(Vi, smooth)
        return Vi

    stat, xe, ye, _ = binned_statistic_2d(x, y, z, statistic="mean", bins=nbin)
    cnt, _, _, _ = binned_statistic_2d(x, y, z, statistic="count", bins=nbin)
    xc = 0.5 * (xe[:-1] + xe[1:]); yc = 0.5 * (ye[:-1] + ye[1:])
    Xc, Yc = np.meshgrid(xc, yc, indexing="ij")
    valid = np.isfinite(stat) & (cnt >= min_count)
    pts = np.column_stack([Xc[valid], Yc[valid]])
    xi = np.linspace(x.min(), x.max(), nx)
    yi = np.linspace(y.min(), y.max(), ny)
    Xi, Yi = np.meshgrid(xi, yi)
    Zi = _fill(stat[valid], pts, Xi, Yi)
    if c is None:
        return Xi, Yi, Zi
    cstat, _, _, _ = binned_statistic_2d(x, y, c, statistic="mean", bins=nbin)
    Ci = _fill(cstat[valid], pts, Xi, Yi)
    return Xi, Yi, Zi, Ci


# ─────────────────────────────────────────────────────────────────────────────
# Plot styling — Times-New-Roman-compatible serif, everything bold
# ─────────────────────────────────────────────────────────────────────────────
def _pick_serif():
    import matplotlib.font_manager as fm
    names = {f.name for f in fm.fontManager.ttflist}
    for want in ("Times New Roman", "Nimbus Roman", "Liberation Serif",
                 "DejaVu Serif"):
        if want in names:
            return want
    return "serif"


SERIF = _pick_serif()


# ─────────────────────────────────────────────────────────────────────────────
# Bright primary colour palette (replaces the earlier muted / dark scheme).
# Drives the default cycle and is referenced by the semantic constants below.
# ─────────────────────────────────────────────────────────────────────────────
BLUE = "#1f77ff"
ORANGE = "#ff8c00"
GREEN = "#00b84d"
RED = "#ee2233"
PURPLE = "#9b30ff"
CYAN = "#12c4e6"
YELLOW = "#ffcc00"
PINK = "#ff3d8b"
TEAL = "#00b3a4"
LIME = "#8ce65a"
AMBER = "#ffb200"
BRIGHT = [BLUE, ORANGE, GREEN, RED, PURPLE, CYAN, YELLOW, PINK, TEAL, LIME]

# Stable colours for the three deterministic confidence components, reused
# across every component-relationship figure (theme 09). C1 = blue, C2 = orange,
# C4 = purple (green/red are reserved for informative/uninformative semantics).
COMPONENT_COLOR = {
    "c1_score": BLUE, "c2_score": ORANGE, "c4_score": PURPLE,
    "mean_c1": BLUE, "mean_c2": ORANGE, "mean_c4": PURPLE,
}

# 21 maximally-distinct bright colours giving every organism a stable identity.
# Deterministic by sorted organism name (see organism_color_map) so the same
# organism keeps the same colour across every figure. Kept vivid so each point
# is separable on white at 300 dpi; a few earthy hues are unavoidable at 21.
ORG_PALETTE = [
    "#e6194b",  # red
    "#f58231",  # orange
    "#d9b600",  # gold
    "#3cb44b",  # green
    "#12b886",  # emerald
    "#00b8d4",  # cyan
    "#1f77ff",  # blue
    "#2b3fd4",  # indigo
    "#911eb4",  # purple
    "#b088ff",  # lavender
    "#f032e6",  # magenta
    "#ff1493",  # deep pink
    "#f58ab0",  # rose
    "#9a6324",  # brown
    "#ff9e4a",  # apricot
    "#94c400",  # yellow-green
    "#469990",  # teal
    "#808000",  # olive
    "#b30000",  # maroon
    "#6e6e6e",  # grey
    "#1a1a1a",  # near-black
]


def apply_style():
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": [SERIF],
        "font.weight": "bold",
        "axes.titleweight": "bold",
        "axes.labelweight": "bold",
        "axes.titlesize": 15,
        "axes.labelsize": 13,
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
        "legend.fontsize": 11,
        "figure.titleweight": "bold",
        "figure.titlesize": 16,
        "savefig.dpi": 300,
        "axes.prop_cycle": plt.cycler(color=BRIGHT),
        "axes.grid": True,
        "grid.alpha": 0.3,
    })


def boldticks(ax):
    for lab in list(ax.get_xticklabels()) + list(ax.get_yticklabels()):
        lab.set_fontweight("bold")


# Consistent colours for the two gene classes / composition categories
COL_INFO = GREEN          # green = informative
COL_UNINFO = RED          # red = uninformative / hypothetical
COL_OPERON = BLUE         # blue = operonic
COL_NONOPERON = ORANGE    # orange = non-operonic
COMP_COLORS = {
    "all_informative": GREEN,
    "majority_informative": LIME,
    "equal": YELLOW,
    "majority_hypothetical": ORANGE,
    "all_hypothetical": RED,
}
COMP_ORDER = ["all_informative", "majority_informative", "equal",
              "majority_hypothetical", "all_hypothetical"]


def organism_color_map(organisms):
    """Deterministic {organism: bright colour}, ordered by sorted full name so
    the same organism keeps the same colour in every figure of the suite."""
    orgs = sorted(set(organisms))
    return {o: ORG_PALETTE[i % len(ORG_PALETTE)] for i, o in enumerate(orgs)}


def add_organism_index(fig, color_map, ncol=7, y=0.94, fontsize=9):
    """Draw a global organism colour index OUTSIDE the axes, centred across the
    top of the figure (below the suptitle). The caller must reserve top margin
    (e.g. fig.subplots_adjust(top=...)) so nothing overlaps. Returns the Legend."""
    from matplotlib.lines import Line2D
    orgs = sorted(color_map)
    handles = [Line2D([0], [0], marker="o", linestyle="none", markersize=8,
                      markerfacecolor=color_map[o], markeredgecolor="black",
                      markeredgewidth=0.4, label=short_label(o)) for o in orgs]
    leg = fig.legend(handles=handles, loc="upper center",
                     bbox_to_anchor=(0.5, y), ncol=ncol, frameon=True,
                     fontsize=fontsize, handletextpad=0.3, columnspacing=1.1,
                     borderpad=0.5, borderaxespad=0.0)
    leg.get_frame().set_edgecolor("black")
    for t in leg.get_texts():
        t.set_fontstyle("italic")
        t.set_fontweight("bold")
    return leg


def savefig(fig, path: Path, dpi=300):
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, bbox_inches="tight", dpi=dpi)
    plt.close(fig)
    print(f"[c3fig] wrote {path}", file=sys.stderr)


def write_tsv(df: pd.DataFrame, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)
    print(f"[c3fig] wrote {path}", file=sys.stderr)


# ─────────────────────────────────────────────────────────────────────────────
# Per-figure runner — keeps each figure script tiny.
# Each figure script defines make(genes, operons, outdir) and calls
# figure_main(make).  Genes/operons come from the pickle caches.
# ─────────────────────────────────────────────────────────────────────────────
import argparse as _argparse


def figure_main(make_fn, subdir=""):
    import inspect
    ap = _argparse.ArgumentParser()
    ap.add_argument("--stats-dir", required=True,
                    help=".../scoring/c3-genes-comprehensive-stats (holds _cache)")
    args = ap.parse_args()
    stats_dir = Path(args.stats_dir)
    cache = stats_dir / "_cache"
    genes = load_cache(cache / "genes.pkl")
    operons = load_cache(cache / "operons.pkl")
    outdir = stats_dir / "figures"
    if subdir:
        outdir = outdir / subdir
    outdir.mkdir(parents=True, exist_ok=True)
    apply_style()
    # Pass adjacent-pairs as a 4th arg only if the figure asks for it.
    nparams = len(inspect.signature(make_fn).parameters)
    if nparams >= 4:
        pairs_path = cache / "adjacent_pairs.pkl"
        pairs = load_cache(pairs_path) if pairs_path.is_file() else build_adjacent_pairs(genes)
        make_fn(genes, operons, outdir, pairs)
    else:
        make_fn(genes, operons, outdir)


# Human genus label for compact organism ticks: "Escherichia coli …" -> "E. coli"
def short_label(organism: str) -> str:
    parts = organism.replace("_", " ").split()
    if len(parts) >= 2:
        return f"{parts[0][0]}. {parts[1]}"
    return organism[:16]


# Compact a long product descriptor for plot labels without losing meaning.
_DESC_SHORTEN = [
    (re.compile(r"\b50S ribosomal protein\b", re.I), "50S-rp"),
    (re.compile(r"\b30S ribosomal protein\b", re.I), "30S-rp"),
    (re.compile(r"\bribosomal protein\b", re.I), "rp"),
    (re.compile(r"\bDNA-directed RNA polymerase\b", re.I), "RNAP"),
    (re.compile(r"\bATP synthase\b", re.I), "ATP-syn"),
    (re.compile(r"\bsubunit\b", re.I), "su"),
    (re.compile(r"\btranscriptional regulator\b", re.I), "transc. reg."),
]


def short_desc(desc: str, maxlen: int = 42) -> str:
    s = str(desc)
    for rx, repl in _DESC_SHORTEN:
        s = rx.sub(repl, s)
    s = re.sub(r"\s+", " ", s).strip()
    if len(s) > maxlen:
        s = s[: maxlen - 1] + "\u2026"
    return s
