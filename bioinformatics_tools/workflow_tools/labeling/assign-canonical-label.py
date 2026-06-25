#!/usr/bin/env python3
"""assign-canonical-label.py — margie_sb phase10 (labeling).

Reads consolidated-merged-all-columns.tsv (the full merge, not any
filtered view -- this needs InterPro sub-database columns and
score/threshold columns that the filtered views deliberately strip out)
and assigns each gene one canonical_label, by walking a fixed trust
hierarchy and taking the first source that clears its quality bar with
an informative description.

Bare script, no container -- pure stdlib, same reasoning as every other
margie_sb consolidation/labeling script.

PRIORITY ORDER (highest confidence first, decision-making tools):
  PGAP > TIGRFAM > HAMAP > NCBIFAM > PIRSF > UNIPROT > PFAM > CDD > KEGG
  > EGGNOG > COG > RAST

  - PGAP/TIGRFAM/HAMAP/NCBIFAM: expert-curated, prokaryote-built-from-
    day-one, equivalog/family-rule-level specificity.
  - PIRSF: whole-protein (not domain-only) HMM classification, ranked
    above UniProt because an HMM generalises better across a divergent
    family member than a single best-alignment hit.
  - UNIPROT: Swiss-Prot's curation is gold-standard, but reached via
    BLAST/DIAMOND best-hit -- gated at >=40% identity.
  - PFAM/CDD: domain-level curation spanning all of life, not
    prokaryote-specific.
  - KEGG: a real, per-KO adaptive score threshold, already enforced
    upstream (merge-all-columns.py's TOOL_ROW_FILTERS).
  - EGGNOG/COG: orthology *inference*, broader/coarser calls.
  - RAST: always-available fallback.

  PANTHER excluded -- eukaryote-curation-skewed. MEROPS/TCDB/DBCAN/
  GeneProp are NOT part of this decision chain -- specialist databases
  that only apply to a subset of genes. They're still computed and kept
  as CONFIRMATORY reference columns (see CONFIRMATORY TOOLS below):
  when one of these narrow, highly specific classifications agrees with
  canonical_label, that's useful corroborating evidence; when it doesn't,
  that's worth a reviewer's (or an LLM's) attention, even though it
  shouldn't drive the primary decision itself.

  EC-number cross-validation is intentionally out of scope for this
  hierarchical decision -- handled as an independent, parallel signal by
  add-ec-consensus.py instead (see that module).

GATES: PFAM/TIGRFAM/PGAP/HAMAP/NCBIfam/PIRSF/CDD/KEGG need no additional
gate -- their search already applied a curated threshold upstream
(--cut_ga/--cut_tc for the HMMER-based ones, tracked in
*_threshold_considered columns; InterProScan's and KofamScan's own
internal cutoffs for the rest). UniProt is the one tool genuinely
ungated upstream -- gated here at >=40% identity. COG/EGGNOG/RAST:
accepted on presence.

MULTI-HIT TOOLS AND MULTI-DOMAIN PROTEINS: PGAP/TIGRFAM/HAMAP/NCBIFAM/
PIRSF/PFAM/CDD/KEGG/EGGNOG/COG/MEROPS/TCDB can have multiple hits per
gene in the merged table ("id1: val1; id2: val2") -- and unlike
"competing predictions for the same role", these are often genuinely
DIFFERENT domains within one multi-domain protein (e.g. an N-terminal
and a C-terminal Pfam domain). So for every one of these tools, BOTH are
kept: the full multi-hit string (every domain, unmodified, in its own
reference column) AND a separate best-hit selection (used for the
labeling decision itself, in its own _best_hit_* column).

Best-hit selection prefers the lowest-e-value hit that is ALSO
informative, not simply the lowest e-value, full stop. A multi-domain
protein can have its literal lowest-e-value hit land on an uninformative
domain (e.g. a "Domain of unknown function") while a different,
slightly-higher-e-value hit on the same protein is a real, named domain.
Picking the literal best-by-e-value in that case would wrongly skip the
tool as "uninformative" and fall through to a lower-priority one, even
though it actually has a perfectly good answer available. The fix: rank
all of a tool's own hits by e-value, but prefer the best-ranked one
that's informative; only fall back to the literal best-by-e-value (even
if uninformative) when NONE of that tool's hits are informative. When
this override happens, gate_note records it explicitly (which hit was
used, which one would have won on e-value alone, and why).
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


def is_uninformative(val: str) -> bool:
    v = val.strip().lower()
    if not v or v in _UNINFORMATIVE:
        return True
    for prefix in _UNINFORMATIVE_PREFIXES:
        if v.startswith(prefix):
            return True
    if _UNKNOWN_FUNC_RE.match(v):
        return True
    if _DB_ID_HYPOTHETICAL_RE.match(v):
        return True
    return False


# ─── STEP 0b: e-value formatting ────────────────────────────────────────────────

def safe_float(s: str, default: float) -> float:
    try:
        return float(s)
    except (TypeError, ValueError):
        return default


def normalize_evalue(value: str) -> str:
    """"4.96051E-66" -> "4.96e-66", "1.300e-55" -> "1.3e-55". Genuine zero
    (HMMER's own underflow result for an extremely significant hit) is
    left exactly as reported -- 0 means "more significant than floating
    point can represent," not missing data."""
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
    """Parse 'id1: val1; id2: val2' into {id1: val1, id2: val2}."""
    result: dict[str, str] = {}
    if not value:
        return result
    for part in value.split("; "):
        if ": " in part:
            k, v = part.split(": ", 1)
            result[k] = v
    return result


class Hit(NamedTuple):
    hit_id: str
    description: str
    stat: str


def parse_all_hits(row: dict[str, str], id_col: str, desc_col: str, stat_col: str) -> list[Hit]:
    """Every hit for this tool, not just the best one. Single-hit columns
    (bare, unwrapped values) return a 1-element list."""
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
    """From all of one tool's hits, select the representative used for
    labeling: the best-ranked (lowest e-value, or highest identity/score
    for rank_mode="max") hit that is ALSO informative, preferred over the
    literal best-ranked hit if that one is uninformative. Falls back to
    the literal best-ranked hit (even if uninformative) only when none of
    this tool's hits are informative.

    Returns (chosen_hit, overrode_literal_best) -- overrode_literal_best
    is True exactly when an uninformative top-ranked hit was passed over
    in favour of a lower-ranked but informative one (e.g. a domain of
    unknown function passed over for a named domain on the same protein).
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
    return ToolEvaluation(tool_name, "", "", False, True, "no hit", False,
                          confirmatory=confirmatory)


# ─── Generic accession-keyed multi-hit evaluator ────────────────────────────────
# Shared by every tool whose hits are accession-keyed and ranked by
# e-value (lower is better): PGAP, TIGRFAM, HAMAP, NCBIFAM, PIRSF, PFAM,
# CDD, KEGG, EGGNOG, COG, MEROPS, TCDB. Each calls this with its own
# column names and an optional extra gate.
def evaluate_accession_tool(
    row: dict[str, str],
    tool_name: str,
    id_col: str,
    desc_col: str,
    stat_col: str,
    gate_note_suffix: str = "",
    confirmatory: bool = False,
) -> ToolEvaluation:
    hits = parse_all_hits(row, id_col, desc_col, stat_col)
    if not hits:
        return _no_hit(tool_name, confirmatory)

    chosen, overrode = select_best_hit(hits, rank_mode="min")
    informative = not is_uninformative(chosen.description)
    stat_norm = normalize_evalue(chosen.stat)

    all_ids = ";".join(h.hit_id for h in hits)
    all_descs = "; ".join(f"{h.hit_id}: {h.description}" for h in hits) if len(hits) > 1 else chosen.description
    all_stats = "; ".join(f"{h.hit_id}: {normalize_evalue(h.stat)}" for h in hits) if len(hits) > 1 else stat_norm

    note = f"evalue={stat_norm}"
    if overrode:
        literal_best = min(hits, key=lambda h: safe_float(h.stat, float("inf")))
        note += (f" (used {chosen.hit_id} over lower-evalue-but-uninformative "
                f"{literal_best.hit_id}, evalue={normalize_evalue(literal_best.stat)})")
    if gate_note_suffix:
        note += f" {gate_note_suffix}"

    return ToolEvaluation(tool_name, chosen.hit_id, chosen.description, informative,
                          True, note, informative, stat_norm,
                          all_ids, all_descs, all_stats, confirmatory)


# ─── STEP 1: PGAP ────────────────────────────────────────────────────────────────
# No gate needed -- PGAP's HMM search already applied its own curated
# trusted-score cutoff at search time (tracked in PGAP_threshold_considered).
def eval_pgap(row: dict[str, str]) -> ToolEvaluation:
    return evaluate_accession_tool(row, "PGAP", "PGAP_id", "PGAP_description", "PGAP_full_seq_evalue")


# ─── STEP 2: TIGRFAM ─────────────────────────────────────────────────────────────
# No gate needed -- same reasoning as PGAP (--cut_tc applied at search time).
def eval_tigrfam(row: dict[str, str]) -> ToolEvaluation:
    return evaluate_accession_tool(row, "TIGRFAM", "TIGRFAM_id", "TIGRFAM_description", "TIGRFAM_full_seq_evalue")


# ─── STEP 3: HAMAP ───────────────────────────────────────────────────────────────
# No gate needed -- InterProScan only emits a HAMAP hit when its own
# internal cutoff is already satisfied.
def eval_hamap(row: dict[str, str]) -> ToolEvaluation:
    return evaluate_accession_tool(row, "HAMAP", "INTERPRO_HAMAP_id", "INTERPRO_HAMAP_description", "INTERPRO_HAMAP_evalue")


# ─── STEP 4: NCBIFAM ─────────────────────────────────────────────────────────────
def eval_ncbifam(row: dict[str, str]) -> ToolEvaluation:
    return evaluate_accession_tool(row, "NCBIFAM", "INTERPRO_NCBIFAM_id", "INTERPRO_NCBIFAM_description", "INTERPRO_NCBIFAM_evalue")


# ─── STEP 5: PIRSF ───────────────────────────────────────────────────────────────
def eval_pirsf(row: dict[str, str]) -> ToolEvaluation:
    return evaluate_accession_tool(row, "PIRSF", "INTERPRO_PIRSF_id", "INTERPRO_PIRSF_description", "INTERPRO_PIRSF_evalue")


# ─── STEP 6: UNIPROT ─────────────────────────────────────────────────────────────
# The one tool genuinely ungated upstream -- explicit gate: identity >= 40%.
# Ranked by HIGHEST identity (not e-value) -- bespoke, not the generic helper.
_UNIPROT_PIDENT_MIN = 40.0


def eval_uniprot(row: dict[str, str]) -> ToolEvaluation:
    hits = parse_all_hits(row, "UNIPROT_id", "UNIPROT_description", "UNIPROT_percent_identity")
    if not hits:
        return _no_hit("UNIPROT")

    # Among hits clearing the identity gate, prefer the informative one
    # with the highest identity; fall back to the single highest-identity
    # hit overall if none qualify.
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
# No gate needed (--cut_ga applied at search time). PFAM hits qualify even
# when uninformative (DUF/UPF) -- a Pfam domain-of-unknown-function hit is
# still a real, specific signal worth keeping as a label of last resort
# within this tool. (select_best_hit already prefers an informative
# domain when this protein has one -- see the module docstring; this
# qualifies=True is only relevant when EVERY one of PFAM's own domains
# for this gene is uninformative.)
def eval_pfam(row: dict[str, str]) -> ToolEvaluation:
    ev = evaluate_accession_tool(row, "PFAM", "PFAM_id", "PFAM_description", "PFAM_full_seq_evalue")
    if ev.hit_id and not ev.qualifies:
        return ev._replace(qualifies=True,
                           gate_note=ev.gate_note + " (uninformative, e.g. DUF -- kept anyway)")
    return ev


# ─── STEP 8: CDD ─────────────────────────────────────────────────────────────────
def eval_cdd(row: dict[str, str]) -> ToolEvaluation:
    return evaluate_accession_tool(row, "CDD", "INTERPRO_CDD_id", "INTERPRO_CDD_description", "INTERPRO_CDD_evalue")


# ─── STEP 9: KEGG ────────────────────────────────────────────────────────────────
# No gate needed here -- merge-all-columns.py's TOOL_ROW_FILTERS already
# dropped every KEGG row with KEGG_is_above_threshold != "true" before this
# script ever sees the data, so any KEGG hit reaching this table already
# cleared its own per-KO adaptive threshold.
def eval_kegg(row: dict[str, str]) -> ToolEvaluation:
    return evaluate_accession_tool(row, "KEGG", "KEGG_id", "KEGG_description", "KEGG_evalue",
                                   gate_note_suffix="(pre-filtered above KofamScan adaptive threshold upstream)")


# ─── STEP 10: EGGNOG ─────────────────────────────────────────────────────────────
# Accepted on presence -- no gate. Now multi-hit-aware (previously this
# tool used a bare row.get(), missing the same informative-preference
# logic every other tool gets).
def eval_eggnog(row: dict[str, str]) -> ToolEvaluation:
    return evaluate_accession_tool(row, "EGGNOG", "EGGNOG_id", "EGGNOG_description", "EGGNOG_evalue")


# ─── STEP 11: COG ────────────────────────────────────────────────────────────────
# Accepted on presence -- no gate. Now multi-hit-aware (same fix as EGGNOG).
def eval_cog(row: dict[str, str]) -> ToolEvaluation:
    return evaluate_accession_tool(row, "COG", "COG_id", "COG_description", "COG_evalue")


# ─── STEP 12: RAST ───────────────────────────────────────────────────────────────
# Always-available fallback -- no id/e-value of its own, accepted on
# presence. Last in priority order: curation consistency varies a lot
# across RAST/SEED's subsystems.
def eval_rast(row: dict[str, str]) -> ToolEvaluation:
    desc = row.get("RAST_description", "")
    if not desc:
        return _no_hit("RAST")
    informative = not is_uninformative(desc)
    return ToolEvaluation("RAST", "", desc, informative, True, "RAST/SEED description",
                          informative)


# Priority order -- the only place the trust hierarchy is actually encoded.
_EVALUATORS: list[Callable[[dict[str, str]], ToolEvaluation]] = [
    eval_pgap, eval_tigrfam, eval_hamap, eval_ncbifam, eval_pirsf,
    eval_uniprot, eval_pfam, eval_cdd, eval_kegg, eval_eggnog, eval_cog, eval_rast,
]


# ─── CONFIRMATORY TOOLS ──────────────────────────────────────────────────────────
# MEROPS/TCDB/dbCAN: narrow, highly specific specialist classifications
# (peptidases, transporters, CAZymes respectively). Never compete for
# canonical_label -- they only apply to a subset of genes, not general
# annotation evidence -- but when one DOES hit, it's valuable corroborating
# (or contradicting) evidence for canonical_label, worth surfacing for a
# reviewer (or an LLM) to cross-check, even though it doesn't drive the
# decision itself.
def eval_merops(row: dict[str, str]) -> ToolEvaluation:
    return evaluate_accession_tool(row, "MEROPS", "MEROPS_id", "MEROPS_description", "MEROPS_evalue",
                                   confirmatory=True)


def eval_tcdb(row: dict[str, str]) -> ToolEvaluation:
    return evaluate_accession_tool(row, "TCDB", "TCDB_id", "TCDB_description", "TCDB_evalue",
                                   confirmatory=True)


def eval_dbcan(row: dict[str, str]) -> ToolEvaluation:
    # dbCAN produces exactly one row per gene (confirmed in its own
    # processing script) -- never multi-hit-wrapped, so no parse_all_hits
    # needed. It also has no e-value of its own; its 3-tool voting
    # consensus count (DBCAN_number_of_tools_hit) is the closest
    # confidence proxy it has.
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


# ─── STEP 13: build the hierarchy + audit trail ─────────────────────────────────
def build_audit_trail(evaluations: list[ToolEvaluation], winner: Optional[ToolEvaluation]) -> str:
    """The exhaustive record: every decision-tool, hit or not, each tagged
    with its outcome and gate note. Answers "why wasn't tool X used" for
    any tool, not just the winner. Confirmatory tools (MEROPS/TCDB/DBCAN)
    are reported separately -- see build_confirmatory_summary()."""
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
    """The full trust-order chain: ALL 12 decision tools, always listed in
    priority order, each marked hit/no-hit -- with the winner's reason
    attached in parentheses at the end."""
    chain = " > ".join(
        f"{ev.tool_name}({'hit' if ev.description else 'no hit'})" for ev in evaluations
    )
    if winner is None:
        return f"{chain} (no qualifying winner -- every hit failed its gate or was uninformative)"
    return f"{chain} (winner: {winner.tool_name}, {winner.gate_note})"


def build_confirmatory_summary(confirmatory_evals: list[ToolEvaluation], canonical_label: str) -> str:
    """For MEROPS/TCDB/DBCAN: report whether each hit, agrees/disagrees/
    is silent relative to canonical_label -- a quick cross-check signal,
    not a decision input. "Agreement" here is a crude substring check
    (shared informative word), good enough to flag for a reviewer/LLM to
    actually judge, not meant to be a final word itself."""
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
    """Returns (canonical_label, label_source, label_source_id,
    label_hierarchy, label_audit_trail, evaluations, confirmatory_evals,
    confirmatory_summary).

    Every decision tool is evaluated regardless of earlier results (no
    short-circuiting) -- the full set of per-tool outcomes is what makes
    this auditable from the output table alone.
    """
    evaluations = [fn(row) for fn in _EVALUATORS]
    confirmatory_evals = [fn(row) for fn in _CONFIRMATORY_EVALUATORS]

    # Pass 1: first evaluation that qualifies (informative + gate passed).
    winner = next((ev for ev in evaluations if ev.qualifies), None)

    # Pass 2: raw fallback -- first tool with ANY description at all,
    # informativeness and gates no longer considered.
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
#
# Output column order:
#   1. organism_name, feature_id, na_seq, aa_seq, na_length, aa_length
#   2. ENVELOPE_envelope_type / ENVELOPE_inference_basis -- the genome-
#      level diderm/monoderm call and, critically, WHY -- e.g. a real
#      marker-based call vs. the conservative no-evidence default used
#      for wall-less edge cases like Mycoplasma (see
#      enrich_with_envelope.py / infer_envelope_type.py).
#   3. domain, gene_id, gene_start, gene_end, RAST_feature_type, RAST_strand
#   4. canonical_label + decision metadata (label_source/_id/_hierarchy/
#      _audit_trail) -- the answer and its full reasoning, grouped together.
#   5. every decision tool's evidence, in priority order: {TOOL}_all_ids/
#      _all_descriptions/_all_stats (every domain/hit, unmodified --
#      reference for multi-domain proteins) followed by {TOOL}_best_hit_id/
#      _best_hit_description/_best_hit_stat (the one actually used, if any).
#   6. RAST_description (RAST is the lowest-priority decision tool).
#   7. confirmatory tools (MEROPS/TCDB/DBCAN): same evidence shape, plus
#      label_confirmatory_summary noting agreement/disagreement with
#      canonical_label.

_IDENTITY_COLUMNS = [
    "organism_name", "feature_id", "na_seq", "aa_seq", "na_length", "aa_length",
    "ENVELOPE_envelope_type", "ENVELOPE_inference_basis",
    "domain", "gene_id", "gene_start", "gene_end",
    "RAST_feature_type", "RAST_strand",
]

# tool_name (matching _EVALUATORS order) -> stat column's real name, used
# only to label the _best_hit_stat column meaningfully (e.g. "evalue").
_STAT_LABELS: dict[str, str] = {
    "PGAP": "evalue", "TIGRFAM": "evalue", "HAMAP": "evalue", "NCBIFAM": "evalue",
    "PIRSF": "evalue", "UNIPROT": "percent_identity", "PFAM": "evalue", "CDD": "evalue",
    "KEGG": "evalue", "EGGNOG": "evalue", "COG": "evalue",
    "MEROPS": "evalue", "TCDB": "evalue", "DBCAN": "tools_agreeing",
}


def _evidence_columns(tool_name: str) -> list[str]:
    stat_label = _STAT_LABELS[tool_name]
    return [
        f"{tool_name}_all_ids", f"{tool_name}_all_descriptions", f"{tool_name}_all_{stat_label}s",
        f"{tool_name}_best_hit_id", f"{tool_name}_best_hit_description", f"{tool_name}_best_hit_{stat_label}",
    ]


_DECISION_TOOL_NAMES = ["PGAP", "TIGRFAM", "HAMAP", "NCBIFAM", "PIRSF", "UNIPROT",
                        "PFAM", "CDD", "KEGG", "EGGNOG", "COG"]
_CONFIRMATORY_TOOL_NAMES = ["MEROPS", "TCDB", "DBCAN"]

_OUTPUT_COLUMNS = (
    _IDENTITY_COLUMNS
    + ["canonical_label", "label_source", "label_source_id", "label_hierarchy", "label_audit_trail"]
    + [col for name in _DECISION_TOOL_NAMES for col in _evidence_columns(name)]
    + ["RAST_description"]
    + [col for name in _CONFIRMATORY_TOOL_NAMES for col in _evidence_columns(name)]
    + ["label_confirmatory_summary"]
)


def _write_evidence(out_row: dict[str, str], ev: ToolEvaluation) -> None:
    stat_label = _STAT_LABELS[ev.tool_name]
    out_row[f"{ev.tool_name}_all_ids"] = ev.all_ids
    out_row[f"{ev.tool_name}_all_descriptions"] = ev.all_descriptions
    out_row[f"{ev.tool_name}_all_{stat_label}s"] = ev.all_stats
    out_row[f"{ev.tool_name}_best_hit_id"] = ev.hit_id
    out_row[f"{ev.tool_name}_best_hit_description"] = ev.description
    out_row[f"{ev.tool_name}_best_hit_{stat_label}"] = ev.stat_value


def main() -> None:
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
            out_row["canonical_label"] = label
            out_row["label_source"] = source
            out_row["label_source_id"] = source_id
            out_row["label_hierarchy"] = hierarchy
            out_row["label_audit_trail"] = trail
            for ev in evaluations:
                if ev.tool_name == "RAST":
                    continue  # no id/stat of its own -- RAST_description alone, written below
                _write_evidence(out_row, ev)
            out_row["RAST_description"] = row.get("RAST_description", "")
            for ev in confirmatory_evals:
                _write_evidence(out_row, ev)
            out_row["label_confirmatory_summary"] = confirmatory_summary

            writer.writerow(out_row)
            by_source[source] = by_source.get(source, 0) + 1
            n += 1

    print(f"[assign-canonical-label] Labeled {n} genes → {output_path}")
    for source, count in sorted(by_source.items(), key=lambda kv: -kv[1]):
        print(f"    {source:10s} {count:6d} ({100.0 * count / n:.1f}%)")


if __name__ == "__main__":
    main()
