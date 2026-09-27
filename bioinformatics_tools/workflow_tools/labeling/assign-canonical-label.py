#!/usr/bin/env python3
"""assign-canonical-label.py — labeling stage: best consensus product descriptor.

Reads consolidated-merged-all-columns.tsv and gives each gene one label by
walking a fixed trust hierarchy and taking the first tool whose best hit is
informative and passes its gate:
  Tier 1: PGAP > NCBIFAM > TIGRFAM   -- curated prokaryotic family rules
  Tier 2: HAMAP > PIRSF > UNIPROT    -- family rule / whole-protein HMM / Swiss-Prot hit
  Tier 3: KEGG | Tier 4: EGGNOG | Tier 5: RAST | Tier 6: PFAM | Tier 7: CDD | Tier 8: COG
  Tier 9: fallback -- the first tool with any description, else "No DB hits"
Only UniProt is gated here (identity >= 40%); the other tools were thresholded
upstream. For multi-hit tools the best hit is the best-ranked informative one,
falling back to the literal best; every hit stays in the _all_* columns.
MEROPS/TCDB/dbCAN are confirmatory only. PANTHER is not used.
"""
from __future__ import annotations
import argparse
import csv
import re
import sys
from pathlib import Path
from typing import Callable, NamedTuple, Optional

csv.field_size_limit(10_000_000)


# ─── STEP 0: uninformative-description filtering ───────────────────────────────

_UNINFORMATIVE = frozenset({
    "", "-", ".", "na", "n/a", "none", "null",
    "unknown", "uncharacterized", "uncharacterised", "putative", "predicted",
    "hypothetical protein",
    "conserved hypothetical protein",
    "conserved protein",
    "conserved domain protein",
    "predicted protein",
    "function unknown",
    "domain of unknown function",
    "general function prediction only",
    "poorly characterized",
    "open reading frame",
})
_NULLISH_UNINFORMATIVE = frozenset({"", "-", ".", "na", "n/a", "none", "null"})

_UNINFORMATIVE_PREFIXES = (
    "domain of unknown function",
    "protein of unknown function",
    "family of unknown function",
    "region of unknown function",
    "repeat of unknown function",
    "module of unknown function",
    "duf",
    "upf",
    "uncharacteri",
    "putative uncharacteri",
    "conserved hypothetical",
    "hypothetical",
    "unknown protein",
    "unknown function",
    "orf",
    "pfam uncharacteri",
)

_UNKNOWN_FUNC_RE = re.compile(
    r'^\s*(?:(?:bacterial|viral|archaeal|eukaryotic|fungal|plant|marine|'
    r'transmembrane|integral membrane|membrane)\s+)?'
    r'(?:domain|protein|family|repeat|region|module)\s+of\s+unknown\s+function',
    re.IGNORECASE,
)

_DB_ID_HYPOTHETICAL_RE = re.compile(
    r'^(?:fig\d+|tigr\d+)[:\s].*(?:hypothetical|conserved hypothetical)',
    re.IGNORECASE,
)

# Further uninformative patterns, kept identical to scoring/score-c1-tool-coverage.py:
# eggNOG UPF families, bare DUF tags, "protein containing domains DUF", and
# mid-string "of unknown function"/"uncharacterized" unless a function word appears.
_UPF_ONLY_RE = re.compile(r'^\s*belongs to the upf\d+', re.IGNORECASE)
# eggNOG's "Psort location ..." descriptions give localisation only, not a product.
_PSORT_ONLY_RE = re.compile(r'^\s*psort location\b', re.IGNORECASE)
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
    """Returns True when a description names no function, using the word lists and regexes above."""
    v = val.strip().lower()
    if not v or v in _UNINFORMATIVE:
        return True
    for prefix in _UNINFORMATIVE_PREFIXES:
        if v.startswith(prefix):
            return True
    if (_UNKNOWN_FUNC_RE.match(v) or _DB_ID_HYPOTHETICAL_RE.match(v)
            or _UPF_ONLY_RE.match(v) or _PSORT_ONLY_RE.match(v) or _BARE_DUF_RE.match(v)
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


# ─── STEP 0b: e-value formatting ────────────────────────────────────────────────

def safe_float(s: str, default: float) -> float:
    """Returns float(s), or default when it cannot be parsed."""
    try:
        return float(s)
    except (TypeError, ValueError):
        return default


def normalize_evalue(value: str) -> str:
    """Formats an e-value to three significant digits ("4.96051E-66" -> "4.96e-66"); zero stays as reported."""
    if not value:
        return value
    try:
        val = float(value)
    except (TypeError, ValueError):
        return value
    if val == 0:
        return value
    mantissa, _, exponent = f"{val:.2e}".partition("e")
    mantissa = mantissa.rstrip("0").rstrip(".")
    exp_sign, exp_digits = exponent[0], exponent[1:].lstrip("0") or "0"
    return f"{mantissa}e{exp_sign}{exp_digits}"


# ─── STEP 0c: multi-hit parsing ─────────────────────────────────────────────────

def parse_keyed(value: str) -> dict[str, str]:
    """Parses 'id1: val1; id2: val2' into {id1: val1, id2: val2}."""
    result: dict[str, str] = {}
    if not value:
        return result
    for part in value.split("; "):
        if ": " in part:
            k, v = part.split(": ", 1)
            result[k] = v
    return result


class Hit(NamedTuple):
    """One tool hit: accession, description and ranking statistic."""
    hit_id: str
    description: str
    stat: str


def parse_all_hits(row: dict[str, str], id_col: str, desc_col: str, stat_col: str) -> list[Hit]:
    """Returns every hit of one tool from its id, description and stat columns."""
    id_val = row.get(id_col, "")
    if not id_val:
        return []
    if ";" not in id_val:
        return [Hit(id_val, row.get(desc_col, ""), row.get(stat_col, ""))]
    ids = id_val.split(";")
    desc_map = parse_keyed(row.get(desc_col, ""))
    stat_map = parse_keyed(row.get(stat_col, ""))
    return [Hit(i, desc_map.get(i, ""), stat_map.get(i, "")) for i in ids]


def select_best_hit(
    hits: list[Hit], rank_mode: str = "min",
) -> tuple[Optional[Hit], bool]:
    """Selects the best-ranked informative hit, falling back to the literal best-ranked hit.

    Ranks by lowest stat (rank_mode="min") or highest ("max"). Returns
    (chosen_hit, overrode_literal_best), True when an uninformative top hit was skipped.
    """
    if not hits:
        return None, False
    default = float("inf") if rank_mode == "min" else float("-inf")
    sign = 1 if rank_mode == "min" else -1
    ranked = sorted(hits, key=lambda h: sign * safe_float(h.stat, default))
    literal_best = ranked[0]
    if not is_uninformative(literal_best.description):
        return literal_best, False
    informative = [h for h in ranked if not is_uninformative(h.description)]
    if informative:
        return informative[0], True
    return literal_best, False


# ─── Evaluation record ──────────────────────────────────────────────────────────

class ToolEvaluation(NamedTuple):
    """One tool's outcome for a gene: chosen hit, gate result and all-hit evidence."""
    tool_name: str
    hit_id: str             # chosen (best-informative) id, used for labeling
    description: str        # chosen description, used for labeling
    informative: bool
    gate_passed: bool
    gate_note: str
    qualifies: bool
    stat_value: str = ""     # chosen hit's stat (e-value/identity), normalised
    all_ids: str = ""        # every hit's id, ";"-joined, unmodified from the merge table
    all_descriptions: str = ""   # every hit's description, raw "id: desc; ..." form
    all_stats: str = ""      # every hit's stat, raw "id: val; ..." form
    confirmatory: bool = False   # True for MEROPS/TCDB/DBCAN -- reference only, never a winner


def _no_hit(tool_name: str, confirmatory: bool = False) -> ToolEvaluation:
    """Returns the evaluation for a tool without hits."""
    return ToolEvaluation(tool_name, "", "", False, True, "no hit", False,
                          confirmatory=confirmatory)


# ─── Generic accession-keyed multi-hit evaluator ────────────────────────────────
def evaluate_accession_tool(
    row: dict[str, str],
    tool_name: str,
    id_col: str,
    desc_col: str,
    stat_col: str,
    gate_note_suffix: str = "",
    confirmatory: bool = False,
) -> ToolEvaluation:
    """Evaluates an accession-keyed tool ranked by e-value: picks the best hit and builds its gate note."""
    hits = parse_all_hits(row, id_col, desc_col, stat_col)
    if not hits:
        return _no_hit(tool_name, confirmatory)

    chosen, overrode = select_best_hit(hits, rank_mode="min")
    stat_norm = normalize_evalue(chosen.stat)

    # Uninformative hits do not qualify in pass 1; a missing description falls
    # back to the accession for pass 2.
    desc_norm = chosen.description.strip().lower()
    if desc_norm and desc_norm not in _NULLISH_UNINFORMATIVE:
        effective_desc = chosen.description
        informative = not is_uninformative(chosen.description)
    else:
        effective_desc = chosen.hit_id
        informative = False

    all_ids = ";".join(h.hit_id for h in hits)
    all_descs = "; ".join(f"{h.hit_id}: {h.description}" for h in hits) if len(hits) > 1 else chosen.description
    all_stats = "; ".join(f"{h.hit_id}: {normalize_evalue(h.stat)}" for h in hits) if len(hits) > 1 else stat_norm

    note = f"evalue={stat_norm}"
    if overrode:
        literal_best = min(hits, key=lambda h: safe_float(h.stat, float("inf")))
        note += (f" (used {chosen.hit_id} over lower-evalue-but-uninformative "
                f"{literal_best.hit_id}, evalue={normalize_evalue(literal_best.stat)})")
    if effective_desc != chosen.description:
        note += f" (description missing; accession used as fallback label)"
    if gate_note_suffix:
        note += f" {gate_note_suffix}"

    return ToolEvaluation(tool_name, chosen.hit_id, effective_desc, informative,
                          True, note, informative, stat_norm,
                          all_ids, all_descs, all_stats, confirmatory)


# ─── STEP 1: PGAP ────────────────────────────────────────────────────────────────
def eval_pgap(row: dict[str, str]) -> ToolEvaluation:
    """Evaluates PGAP hits (trusted cutoff applied at search time)."""
    return evaluate_accession_tool(row, "PGAP", "PGAP_id", "PGAP_description", "PGAP_full_seq_evalue")


# ─── STEP 2: TIGRFAM ─────────────────────────────────────────────────────────────
def eval_tigrfam(row: dict[str, str]) -> ToolEvaluation:
    """Evaluates TIGRFAM hits (--cut_tc applied at search time)."""
    return evaluate_accession_tool(row, "TIGRFAM", "TIGRFAM_id", "TIGRFAM_description", "TIGRFAM_full_seq_evalue")


# ─── STEP 3: HAMAP ───────────────────────────────────────────────────────────────
def eval_hamap(row: dict[str, str]) -> ToolEvaluation:
    """Evaluates InterPro HAMAP hits (InterProScan cutoff applied upstream)."""
    return evaluate_accession_tool(row, "HAMAP", "INTERPRO_HAMAP_id", "INTERPRO_HAMAP_description", "INTERPRO_HAMAP_evalue")


# ─── STEP 4: NCBIFAM ─────────────────────────────────────────────────────────────
def eval_ncbifam(row: dict[str, str]) -> ToolEvaluation:
    """Evaluates InterPro NCBIfam hits."""
    return evaluate_accession_tool(row, "NCBIFAM", "INTERPRO_NCBIFAM_id", "INTERPRO_NCBIFAM_description", "INTERPRO_NCBIFAM_evalue")


# ─── STEP 5: PIRSF ───────────────────────────────────────────────────────────────
def eval_pirsf(row: dict[str, str]) -> ToolEvaluation:
    """Evaluates InterPro PIRSF hits."""
    return evaluate_accession_tool(row, "PIRSF", "INTERPRO_PIRSF_id", "INTERPRO_PIRSF_description", "INTERPRO_PIRSF_evalue")


# ─── STEP 6: UNIPROT ─────────────────────────────────────────────────────────────
# Gated here at identity >= 40% and ranked by highest identity.
_UNIPROT_PIDENT_MIN = 40.0


def eval_uniprot(row: dict[str, str]) -> ToolEvaluation:
    """Evaluates UniProt hits: prefers the highest-identity informative hit above the identity gate."""
    hits = parse_all_hits(row, "UNIPROT_id", "UNIPROT_description", "UNIPROT_percent_identity")
    if not hits:
        return _no_hit("UNIPROT")

    ranked = sorted(hits, key=lambda h: -safe_float(h.stat, float("-inf")))
    gated = [h for h in ranked if safe_float(h.stat, -1.0) >= _UNIPROT_PIDENT_MIN]
    pool = gated if gated else ranked
    informative_pool = [h for h in pool if not is_uninformative(h.description)]
    chosen = informative_pool[0] if informative_pool else pool[0]

    pident = safe_float(chosen.stat, -1.0)
    gate_passed = pident < 0 or pident >= _UNIPROT_PIDENT_MIN
    gate_note = (f"pident={pident:.1f}% {'>=' if gate_passed else '<'} {_UNIPROT_PIDENT_MIN:.0f}% required"
                if pident >= 0 else "pident unavailable, gate skipped")
    informative = not is_uninformative(chosen.description)

    all_ids = ";".join(h.hit_id for h in hits)
    all_descs = "; ".join(f"{h.hit_id}: {h.description}" for h in hits) if len(hits) > 1 else chosen.description
    all_stats = "; ".join(f"{h.hit_id}: {h.stat}" for h in hits) if len(hits) > 1 else chosen.stat

    return ToolEvaluation("UNIPROT", chosen.hit_id, chosen.description, informative,
                          gate_passed, gate_note, informative and gate_passed, chosen.stat,
                          all_ids, all_descs, all_stats)


# ─── STEP 7: PFAM ────────────────────────────────────────────────────────────────
def eval_pfam(row: dict[str, str]) -> ToolEvaluation:
    """Evaluates PFAM hits (--cut_ga upstream); a hit qualifies even when all its domains are DUF/UPF."""
    ev = evaluate_accession_tool(row, "PFAM", "PFAM_id", "PFAM_description", "PFAM_full_seq_evalue")
    if ev.hit_id and not ev.qualifies:
        return ev._replace(qualifies=True,
                           gate_note=ev.gate_note + " (uninformative, e.g. DUF -- kept anyway)")
    return ev


# ─── STEP 8: CDD ─────────────────────────────────────────────────────────────────
def eval_cdd(row: dict[str, str]) -> ToolEvaluation:
    """Evaluates InterPro CDD hits."""
    return evaluate_accession_tool(row, "CDD", "INTERPRO_CDD_id", "INTERPRO_CDD_description", "INTERPRO_CDD_evalue")


# ─── STEP 9: KEGG ────────────────────────────────────────────────────────────────
def eval_kegg(row: dict[str, str]) -> ToolEvaluation:
    """Evaluates KEGG hits; merge-all-columns.py already kept only rows above the KofamScan threshold."""
    return evaluate_accession_tool(row, "KEGG", "KEGG_id", "KEGG_description", "KEGG_evalue",
                                   gate_note_suffix="(pre-filtered above KofamScan adaptive threshold upstream)")


# ─── STEP 10: EGGNOG ─────────────────────────────────────────────────────────────
def eval_eggnog(row: dict[str, str]) -> ToolEvaluation:
    """Evaluates eggNOG hits, accepted on presence."""
    return evaluate_accession_tool(row, "EGGNOG", "EGGNOG_id", "EGGNOG_description", "EGGNOG_evalue")


# ─── STEP 11: COG ────────────────────────────────────────────────────────────────
def eval_cog(row: dict[str, str]) -> ToolEvaluation:
    """Evaluates COG hits, accepted on presence."""
    return evaluate_accession_tool(row, "COG", "COG_id", "COG_description", "COG_evalue")


# ─── RAST ────────────────────────────────────────────────────────────────────────
def eval_rast(row: dict[str, str]) -> ToolEvaluation:
    """Evaluates the RAST/SEED description (no id or e-value), accepted on presence."""
    desc = row.get("RAST_description", "")
    if not desc:
        return _no_hit("RAST")
    informative = not is_uninformative(desc)
    return ToolEvaluation("RAST", "", desc, informative, True, "RAST/SEED description",
                          informative)


# Trust hierarchy in priority order (tiers 1-8); tier 9 is assign_label's pass 2.
_EVALUATORS: list[Callable[[dict[str, str]], ToolEvaluation]] = [
    eval_pgap, eval_ncbifam, eval_tigrfam, eval_hamap, eval_pirsf,
    eval_uniprot, eval_kegg, eval_eggnog, eval_rast, eval_pfam, eval_cdd, eval_cog,
]


# ─── CONFIRMATORY TOOLS ──────────────────────────────────────────────────────────
# MEROPS/TCDB/dbCAN (peptidases, transporters, CAZymes) never pick the label;
# their hits are reported for cross-checking.
def eval_merops(row: dict[str, str]) -> ToolEvaluation:
    """Evaluates MEROPS hits as confirmatory evidence."""
    return evaluate_accession_tool(row, "MEROPS", "MEROPS_id", "MEROPS_description", "MEROPS_evalue",
                                   confirmatory=True)


def eval_tcdb(row: dict[str, str]) -> ToolEvaluation:
    """Evaluates TCDB hits as confirmatory evidence."""
    return evaluate_accession_tool(row, "TCDB", "TCDB_id", "TCDB_description", "TCDB_evalue",
                                   confirmatory=True)


def eval_dbcan(row: dict[str, str]) -> ToolEvaluation:
    """Evaluates the single dbCAN row as confirmatory evidence; its stat is the number of agreeing sub-tools."""
    desc = row.get("DBCAN_description", "")
    if not desc:
        return _no_hit("DBCAN", confirmatory=True)
    hid = row.get("DBCAN_id", "")
    n_tools = row.get("DBCAN_number_of_tools_hit", "")
    informative = not is_uninformative(desc)
    return ToolEvaluation("DBCAN", hid, desc, informative, True,
                          f"{n_tools} sub-tools agree", informative, n_tools,
                          hid, desc, n_tools, confirmatory=True)


_CONFIRMATORY_EVALUATORS: list[Callable[[dict[str, str]], ToolEvaluation]] = [
    eval_merops, eval_tcdb, eval_dbcan,
]


_FALLBACK_NO_HITS = "No DB hits"


# ─── Hierarchy and audit trail ──────────────────────────────────────────────────
def build_audit_trail(evaluations: list[ToolEvaluation], winner: Optional[ToolEvaluation]) -> str:
    """Lists every decision tool with its outcome tag and gate note, " | "-joined."""
    trail_parts = []
    for ev in evaluations:
        tag = "WINNER" if winner is ev else (
            "no_hit" if not ev.description else
            "uninformative" if not ev.informative else
            "gate_failed" if not ev.gate_passed else
            "available_but_lower_priority"
        )
        if ev.description:
            trail_parts.append(f"{ev.tool_name}=[{tag}; {ev.gate_note}]")
        else:
            trail_parts.append(f"{ev.tool_name}=[{tag}]")
    return " | ".join(trail_parts)


def build_hierarchy(evaluations: list[ToolEvaluation], winner: Optional[ToolEvaluation]) -> str:
    """Returns the priority-ordered tool chain marked hit/no hit, with the winner and its gate note."""
    chain = " > ".join(
        f"{ev.tool_name}({'hit' if ev.description else 'no hit'})" for ev in evaluations
    )
    if winner is None:
        return f"{chain} (no qualifying winner -- every hit failed its gate or was uninformative)"
    return f"{chain} (winner: {winner.tool_name}, {winner.gate_note})"


def build_confirmatory_summary(confirmatory_evals: list[ToolEvaluation], canonical_label: str) -> str:
    """Summarises MEROPS/TCDB/dbCAN hits against the label by shared 4+-letter words (regex)."""
    parts = []
    label_words = set(re.findall(r"[a-z]{4,}", canonical_label.lower()))
    for ev in confirmatory_evals:
        if not ev.description:
            parts.append(f"{ev.tool_name}=[no hit]")
            continue
        desc_words = set(re.findall(r"[a-z]{4,}", ev.description.lower()))
        overlap = label_words & desc_words
        relation = "possible_agreement" if overlap else "unclear -- needs review"
        parts.append(f"{ev.tool_name}=[hit; {ev.description[:60]}; {relation}]")
    return " | ".join(parts)


def assign_label(
    row: dict[str, str],
) -> tuple[str, str, str, str, str, list[ToolEvaluation], list[ToolEvaluation], str]:
    """Evaluates every tool and picks the label: first qualifying tool, else first with any description.

    Returns (label, source, source_id, hierarchy, audit_trail, evaluations,
    confirmatory_evals, confirmatory_summary).
    """
    evaluations = [fn(row) for fn in _EVALUATORS]
    confirmatory_evals = [fn(row) for fn in _CONFIRMATORY_EVALUATORS]

    # Pass 1: first evaluation that qualifies (informative + gate passed).
    winner = next((ev for ev in evaluations if ev.qualifies), None)

    # Pass 2: first tool with any description, ignoring informativeness and gates.
    if winner is None:
        winner = next((ev for ev in evaluations if ev.description), None)

    label = winner.description[:200].strip() if winner else _FALLBACK_NO_HITS
    source = winner.tool_name if winner else "NONE"
    source_id = winner.hit_id if winner else ""

    hierarchy = build_hierarchy(evaluations, winner)
    trail = build_audit_trail(evaluations, winner)
    confirmatory_summary = build_confirmatory_summary(confirmatory_evals, label)

    return label, source, source_id, hierarchy, trail, evaluations, confirmatory_evals, confirmatory_summary


# ─── CLI ──────────────────────────────────────────────────────────────────────
# Output columns: identity and envelope, the label with its decision metadata,
# each decision tool's all-hit and best-hit evidence, RAST_description, then
# confirmatory tools and their summary.

_IDENTITY_COLUMNS = [
    "organism_name", "feature_id", "na_seq", "aa_seq", "na_length", "aa_length",
    "ENVELOPE_envelope_type", "ENVELOPE_inference_basis",
    "domain", "gene_id", "gene_start", "gene_end",
    "RAST_feature_type", "RAST_strand",
]

# tool_name -> name used for its stat columns (e.g. "evalue").
_STAT_LABELS: dict[str, str] = {
    "PGAP": "evalue", "TIGRFAM": "evalue", "HAMAP": "evalue", "NCBIFAM": "evalue",
    "PIRSF": "evalue", "UNIPROT": "percent_identity", "PFAM": "evalue", "CDD": "evalue",
    "KEGG": "evalue", "EGGNOG": "evalue", "COG": "evalue",
    "MEROPS": "evalue", "TCDB": "evalue", "DBCAN": "tools_agreeing",
}


def _evidence_columns(tool_name: str) -> list[str]:
    """Returns a tool's six evidence column names."""
    stat_label = _STAT_LABELS[tool_name]
    return [
        f"{tool_name}_all_ids", f"{tool_name}_all_descriptions", f"{tool_name}_all_{stat_label}s",
        f"{tool_name}_best_hit_id", f"{tool_name}_best_hit_description", f"{tool_name}_best_hit_{stat_label}",
    ]


_DECISION_TOOL_NAMES = ["PGAP", "NCBIFAM", "TIGRFAM", "HAMAP", "PIRSF", "UNIPROT",
                        "KEGG", "EGGNOG", "PFAM", "CDD", "COG"]
_CONFIRMATORY_TOOL_NAMES = ["MEROPS", "TCDB", "DBCAN"]

_OUTPUT_COLUMNS = (
    _IDENTITY_COLUMNS
    + [
        "best_consensus_product_descriptor",
        "product_descriptor_source",
        "product_descriptor_source_id",
        "product_descriptor_hierarchy",
        "product_descriptor_audit_trail",
    ]
    + [col for name in _DECISION_TOOL_NAMES for col in _evidence_columns(name)]
    + ["RAST_description"]
    + [col for name in _CONFIRMATORY_TOOL_NAMES for col in _evidence_columns(name)]
    + ["product_descriptor_confirmatory_summary"]
)


def _write_evidence(out_row: dict[str, str], ev: ToolEvaluation) -> None:
    """Writes a tool's evidence columns into the output row."""
    stat_label = _STAT_LABELS[ev.tool_name]
    out_row[f"{ev.tool_name}_all_ids"] = ev.all_ids
    out_row[f"{ev.tool_name}_all_descriptions"] = ev.all_descriptions
    out_row[f"{ev.tool_name}_all_{stat_label}s"] = ev.all_stats
    out_row[f"{ev.tool_name}_best_hit_id"] = ev.hit_id
    out_row[f"{ev.tool_name}_best_hit_description"] = ev.description
    out_row[f"{ev.tool_name}_best_hit_{stat_label}"] = ev.stat_value


def main() -> None:
    """Streams the merged table with csv, labels each gene and writes labeled-genes.tsv."""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", required=True,
                        help="merge-all-columns.py's output TSV (consolidated-merged-all-columns.tsv)")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    input_path = Path(args.input)
    if not input_path.is_file():
        print(f"[assign-canonical-label] ERROR: input not found: {input_path}", file=sys.stderr)
        raise SystemExit(1)

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    by_source: dict[str, int] = {}
    n = 0
    with open(input_path, newline="") as fh, open(output_path, "w", newline="") as out_fh:
        reader = csv.DictReader(fh, delimiter="\t")
        writer = csv.DictWriter(out_fh, fieldnames=_OUTPUT_COLUMNS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in reader:
            (label, source, source_id, hierarchy, trail,
             evaluations, confirmatory_evals, confirmatory_summary) = assign_label(row)

            out_row = {col: row.get(col, "") for col in _IDENTITY_COLUMNS}
            out_row["best_consensus_product_descriptor"] = label
            out_row["product_descriptor_source"] = source
            out_row["product_descriptor_source_id"] = source_id
            out_row["product_descriptor_hierarchy"] = hierarchy
            out_row["product_descriptor_audit_trail"] = trail
            for ev in evaluations:
                if ev.tool_name == "RAST":
                    continue  # RAST has only RAST_description, written below
                _write_evidence(out_row, ev)
            out_row["RAST_description"] = row.get("RAST_description", "")
            for ev in confirmatory_evals:
                _write_evidence(out_row, ev)
            out_row["product_descriptor_confirmatory_summary"] = confirmatory_summary

            writer.writerow(out_row)
            by_source[source] = by_source.get(source, 0) + 1
            n += 1

    print(f"[assign-canonical-label] Labeled {n} genes → {output_path}")
    for source, count in sorted(by_source.items(), key=lambda kv: -kv[1]):
        print(f"    {source:10s} {count:6d} ({100.0 * count / n:.1f}%)")


if __name__ == "__main__":
    main()
