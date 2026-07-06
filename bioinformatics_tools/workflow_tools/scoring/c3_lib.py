#!/usr/bin/env python3
"""Descriptor gate and per-organism loader for the C3 operon-context scorer.

Two responsibilities:
  * is_uninformative() / clean_descriptor() -- decide whether a gene carries a
    real function name and normalise that name to its functional identity (strip
    the leading "SOURCE:" tag and the "raw ## human-readable" duplication). The
    uninformative gate mirrors the authoritative definition in
    score-c1-tool-coverage.py; keep the two consistent.
  * load_organism() -- JOIN one organism's labeled-genes.tsv (coordinates,
    descriptor, aa_seq) with its labeled-genes-operon-info.tsv (operon_id,
    member_count, position, probability) on feature_id into one row per gene.
"""
import csv
import hashlib
import re
from pathlib import Path

import pandas as pd

csv.field_size_limit(10_000_000)


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


def genome_fingerprint(aa_hashes) -> str:
    """Content identity for a genome: sha256 over its SORTED per-gene aa_hashes.

    Order-independent and stable across re-annotation (same proteome -> same
    fingerprint), so it identifies a genome by its CONTENT rather than by the
    user-chosen, arbitrary organism/file name.  This is the key used to detect
    whether a candidate genome is already in the OCC reference (for leave-one-out
    scoring) and to dedupe on add.  Empty hashes are ignored; "" if none.

    Pass an iterable of aa_hash strings, e.g. load_organism(..., compute_hash=
    True)["aa_hash"]."""
    hs = sorted(h for h in aa_hashes if h)
    if not hs:
        return ""
    return hashlib.sha256("\n".join(hs).encode("utf-8")).hexdigest()


# columns of a genome's pool-stats record (order used by the sidecar TSV)
POOL_STAT_FIELDS = ("total_genes", "operonic_genes", "singleton_genes",
                    "n_operons", "n_informative_operons", "n_uninformative_operons")


def genome_pool_stats(genes) -> dict:
    """Descriptive pool statistics for ONE genome, from a load_organism() frame.

    An operon is INFORMATIVE (== the OCC 'qualifying' gate in c3_occ) iff its
    informative members strictly outnumber its uninformative ones; otherwise it
    is UNINFORMATIVE ('disinformative' -- dominated by hypotheticals). Returns a
    dict over POOL_STAT_FIELDS: total genes, genes in any operon, non-operonic
    (singleton) genes, and operon counts split informative / uninformative."""
    total = len(genes)
    op = genes[genes["operon_id"].astype(str).str.startswith("operon_")]
    n_op = n_inf_op = n_unf_op = 0
    for _, sub in op.groupby("operon_id"):
        inf = sum(1 for u, d in zip(sub["uninformative"], sub["clean_descriptor"])
                  if (not bool(u)) and str(d).strip())
        n_op += 1
        if inf > (len(sub) - inf) and inf > 0:
            n_inf_op += 1
        else:
            n_unf_op += 1
    return dict(total_genes=total, operonic_genes=len(op),
                singleton_genes=total - len(op), n_operons=n_op,
                n_informative_operons=n_inf_op, n_uninformative_operons=n_unf_op)


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
