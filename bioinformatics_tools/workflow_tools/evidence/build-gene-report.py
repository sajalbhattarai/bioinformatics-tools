#!/usr/bin/env python3
"""build-gene-report.py -- margie_sb phase14 (evidence): one fully tabulated,
self-contained GENE ANNOTATION REPORT per protein-coding gene, read straight
off consolidation's own current column names.

Runs after phase11 (scoring) and phase12/13 (fingerprint/synteny), before
phase15 (llm) -- pure CPU/file-IO, its output is useful and inspectable on
its own with or without ever calling a model.

Format follows the user's own sample-report-before-llm.txt: one shared
11-column table (Tool/ID/Description/Domain/Category/EC/Score/Bitscore/
%Identity/E-value/Other) for every general+specialized database, and one
shared 6-column table (Tool/Prediction/Domain/Score-Probability/
Topology-CleavageSite/Other) for every localization tool -- every
registered tool always gets a row, "-" where a column is genuinely not
sourceable for that tool (real per-tool data backs each column; see
individual row builders for which fields exist and which are real "-").

Two deliberate deviations from a literal reading of that sample, both
because real column-availability disagrees with the hand-typed dash
pattern (the sample was written before checking real data; not every
cell's dash/blank guess was right):
  - NCBIFAM has no Score/Bitscore/%Identity (InterPro's NCBIFAM member-DB
    columns only carry id/description/evalue/match coordinates).
  - DBCAN has no Score/Bitscore/%Identity/evalue at all (confirmed against
    the raw table); its EC column IS real (DBCAN_ec_numbers) even though
    the sample marked it "-".

OPERON MEMBERS is capped to OPERON_WINDOW nearest neighbors on each side of
the candidate (not literally "every remaining operonic gene") -- full
nesting made report size scale with operon size; see this module's own
git history / prior conversation for the 1452-line example that motivated
the cap.
"""
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

csv.field_size_limit(10_000_000)

import pandas as pd

BANNER = "=" * 80
LABEL_W = 24
OPERON_WINDOW = 3  # genes shown upstream and downstream of the candidate
_NO_REAL_OPERON = ("NOT_IN_AN_OPERON", "NOT_APPLICABLE_NON_CODING")


def _clean(val) -> str:
    if val is None:
        return ""
    s = str(val).strip()
    return "" if s.lower() in ("nan", "none", "") else s


def _fmt_pct(val) -> str:
    v = _clean(val)
    if not v:
        return ""
    try:
        return f"{float(v):.1f}%"
    except ValueError:
        return v


def _fmt_eval(val) -> str:
    v = _clean(val)
    if not v:
        return ""
    try:
        return f"{float(v):.2e}"
    except ValueError:
        return v


def _has_value(val) -> bool:
    """True unless empty or the tool's own "-" placeholder (KEGG/EGGNOG/
    DBCAN all use a literal "-" for "checked, nothing found")."""
    return _clean(val) not in ("", "-")


def _int_str(val) -> str:
    """OPERON_gene_position_in_operon is stored as a float-string ("2.0")
    -- display it as a plain integer."""
    v = _clean(val)
    try:
        return str(int(float(v)))
    except ValueError:
        return v


def kv(pairs: list[tuple[str, str]], indent: str = "") -> str:
    return "\n".join(f"{indent}{label:<{LABEL_W}}: {value}" for label, value in pairs)


def pipe_table(header: list[str], rows: list[list[str]], indent: str = "",
               first_col_w: int = 10) -> str:
    """Pipe-delimited table, first column padded for a tidy left edge,
    everything after it plain " | "-joined (matches the user's own sample
    report exactly -- only the leading Tool/Variant/Component column is
    padded there too). Dash separator length matches the rendered header."""
    head_line = indent + " | ".join([header[0].ljust(first_col_w)] + header[1:])
    sep = indent + "-" * (len(head_line) - len(indent))
    lines = [head_line, sep]
    for r in rows:
        cells = [str(r[0]).ljust(first_col_w)] + [str(c) if c else "-" for c in r[1:]]
        lines.append(indent + " | ".join(cells))
    return "\n".join(lines)


# ---- multi-domain field splitting -----------------------------------------
# Consolidation packs >1 domain hit per gene+DB as "ID1: val1; ID2: val2" in
# every column (alignment_from/to, score, evalue, ...) keyed by the same
# domain ID found in the *_id column itself; a single hit is stored bare
# (no "ID: " prefix at all). Both shapes need to come back out as one row
# per domain rather than one semicolon-packed cell, so an LLM (or a human)
# doesn't have to re-split anything itself.
def _parse_id_value_list(s: str) -> list[tuple[str, str]]:
    pairs = []
    for entry in s.split("; "):
        if ": " in entry:
            dom_id, val = entry.split(": ", 1)
            pairs.append((dom_id.strip(), val.strip()))
        else:
            pairs.append(("", entry.strip()))
    return pairs


def split_multi(row: dict, id_col: str, *value_cols: str) -> list[list[str]]:
    """Returns one row per domain: [id, value_col_1, value_col_2, ...].
    Falls back to broadcasting a bare (non-"ID: val") value across every
    domain ID when a column wasn't itself domain-keyed."""
    id_val = _clean(row.get(id_col, ""))
    if not id_val:
        return []
    ids = [i.strip() for i in id_val.split(";") if i.strip()]
    out = [[i] for i in ids]
    for col in value_cols:
        raw = _clean(row.get(col, ""))
        if not raw:
            for r in out:
                r.append("")
            continue
        pairs = _parse_id_value_list(raw)
        if len(pairs) == len(ids):
            for r, (_, v) in zip(out, pairs):
                r.append(v)
        else:
            for r in out:
                r.append(raw)
    return out


# ---- per-tool row builders --------------------------------------------------
# Every function returns rows shaped [Tool, ID, Description, Domain,
# Category, EC, Score, Bitscore, %Identity, E-value, Other] -- one shared
# header for every general/specialized database. A tool that found nothing
# for this gene still returns exactly one all-"-" row (every registered
# tool always appears, never silently omitted).
def _row(tool, id_="", desc="", domain="", category="", ec="", score="",
         bitscore="", pct_identity="", evalue="", other="") -> list[str]:
    return [tool, id_, desc, domain, category, ec, score, bitscore, pct_identity, evalue, other]


def rows_rast(row: dict) -> list[list[str]]:
    desc = _clean(row.get("RAST_description"))
    if not desc:
        return [_row("RAST")]
    subsys = _clean(row.get("RASTTK_subsystem_names"))
    return [_row("RAST", desc=desc, other=f"subsystem: {subsys}" if subsys else "")]


def rows_pgap(row: dict) -> list[list[str]]:
    domains = split_multi(row, "PGAP_id", "PGAP_description", "PGAP_alignment_from",
                           "PGAP_alignment_to", "PGAP_full_seq_evalue")
    if not domains:
        return [_row("PGAP")]
    return [_row("PGAP", id_=i, desc=d, domain=f"{f}-{t}" if f or t else "", evalue=_fmt_eval(e))
            for i, d, f, t, e in domains]


def rows_tigrfam(row: dict) -> list[list[str]]:
    domains = split_multi(row, "TIGRFAM_id", "TIGRFAM_description", "TIGRFAM_alignment_from",
                           "TIGRFAM_alignment_to", "TIGRFAM_full_seq_score", "TIGRFAM_full_seq_evalue")
    if not domains:
        return [_row("TIGRFAM")]
    return [_row("TIGRFAM", id_=i, desc=d, domain=f"{f}-{t}" if f or t else "", score=s, evalue=_fmt_eval(e))
            for i, d, f, t, s, e in domains]


def rows_ncbifam(row: dict) -> list[list[str]]:
    # Via InterPro's NCBIFAM member-DB columns -- no Score/Bitscore/
    # %Identity exist for it at all (only id/description/evalue/match
    # coordinates), unlike the sample's dash pattern, which left Score
    # blank under the assumption it might exist.
    domains = split_multi(row, "INTERPRO_NCBIFAM_id", "INTERPRO_NCBIFAM_description",
                           "INTERPRO_NCBIFAM_match_start", "INTERPRO_NCBIFAM_match_end",
                           "INTERPRO_NCBIFAM_evalue")
    if not domains:
        return [_row("NCBIFAM")]
    return [_row("NCBIFAM", id_=i, desc=d, domain=f"{f}-{t}" if f or t else "", evalue=_fmt_eval(e))
            for i, d, f, t, e in domains]


def rows_cog(row: dict) -> list[list[str]]:
    id_, desc = _clean(row.get("COG_id")), _clean(row.get("COG_description"))
    if not id_:
        return [_row("COG")]
    return [_row("COG", id_=id_, desc=desc, category=_clean(row.get("COG_func_letter")),
                 pct_identity=_fmt_pct(row.get("COG_identity")), evalue=_fmt_eval(row.get("COG_evalue")))]


def rows_pfam(row: dict) -> list[list[str]]:
    domains = split_multi(row, "PFAM_id", "PFAM_description", "PFAM_alignment_from",
                           "PFAM_alignment_to", "PFAM_full_seq_score", "PFAM_full_seq_evalue")
    if not domains:
        return [_row("PFAM")]
    return [_row("PFAM", id_=i, desc=d, domain=f"{f}-{t}" if f or t else "", score=s, evalue=_fmt_eval(e))
            for i, d, f, t, s, e in domains]


def rows_geneprop(row: dict) -> list[list[str]]:
    domains = split_multi(row, "GENEPROP_id", "GENEPROP_description", "GENEPROP_status")
    if not domains:
        return [_row("GENEPROP")]
    return [_row("GENEPROP", id_=i, desc=d, other=f"status: {s}" if s else "") for i, d, s in domains]


def rows_interpro(row: dict) -> list[list[str]]:
    # No per-accession domain-boundary column exists at this merged,
    # deduplicated level (lives only in the 20 INTERPRO_<memberdb>_*
    # column groups, keyed differently per member DB). GO terms are
    # deduplicated IDs only; INTERPRO_pathways is almost entirely
    # human/animal Reactome IDs inherited from cross-species InterPro
    # mappings -- not meaningful evidence for a bacterial/archaeal gene,
    # dropped entirely rather than dumping hundreds of irrelevant
    # cross-species pathway IDs into every report.
    domains = split_multi(row, "INTERPRO_id", "INTERPRO_description")
    if not domains:
        return [_row("INTERPRO")]
    go_ids = sorted(set(re.findall(r"GO:\d+", _clean(row.get("INTERPRO_go_terms")))))
    go_stat = f"GO: {','.join(go_ids)}" if go_ids else ""
    out = [_row("INTERPRO", id_=i, desc=d) for i, d in domains]
    if go_stat:
        out[0][-1] = go_stat
    return out


_EC_RE = re.compile(r"\[EC:([^\]]+)\]")


def rows_kegg(row: dict) -> list[list[str]]:
    id_, desc = _clean(row.get("KEGG_id")), _clean(row.get("KEGG_description"))
    if not id_:
        return [_row("KEGG")]
    ec_match = _EC_RE.search(desc)
    other_parts = []
    if _has_value(row.get("KEGG_reaction_ids")):
        other_parts.append(f"reactions: {_clean(row.get('KEGG_reaction_ids'))}")
    if _has_value(row.get("KEGG_compound_ids")):
        other_parts.append(f"compounds: {_clean(row.get('KEGG_compound_ids'))}")
    return [_row("KEGG", id_=id_, desc=desc, ec=ec_match.group(1) if ec_match else "",
                 evalue=_fmt_eval(row.get("KEGG_evalue")), other=", ".join(other_parts))]


def rows_eggnog(row: dict) -> list[list[str]]:
    id_, desc = _clean(row.get("EGGNOG_id")), _clean(row.get("EGGNOG_description"))
    if not id_:
        return [_row("EGGNOG")]
    other_parts = []
    if _has_value(row.get("EGGNOG_GO_terms")):
        other_parts.append(f"GO: {_clean(row.get('EGGNOG_GO_terms'))}")
    if _has_value(row.get("EGGNOG_KEGG_ko")):
        other_parts.append(f"KEGG_ko: {_clean(row.get('EGGNOG_KEGG_ko'))}")
    if _has_value(row.get("EGGNOG_PFAMs")):
        other_parts.append(f"PFAMs: {_clean(row.get('EGGNOG_PFAMs'))}")
    return [_row("EGGNOG", id_=id_, desc=desc, category=_clean(row.get("EGGNOG_COG_category")),
                 evalue=_fmt_eval(row.get("EGGNOG_evalue")), other=", ".join(other_parts))]


def rows_uniprot(row: dict) -> list[list[str]]:
    id_, desc = _clean(row.get("UNIPROT_id")), _clean(row.get("UNIPROT_description"))
    if not id_:
        return [_row("UNIPROT")]
    other_parts = []
    if _clean(row.get("UNIPROT_gene_name")):
        other_parts.append(f"gene: {_clean(row.get('UNIPROT_gene_name'))}")
    if _clean(row.get("UNIPROT_source_organism")):
        other_parts.append(f"organism: {_clean(row.get('UNIPROT_source_organism'))}")
    return [_row("UNIPROT", id_=id_, desc=desc, bitscore=_clean(row.get("UNIPROT_bitscore")),
                 pct_identity=_fmt_pct(row.get("UNIPROT_percent_identity")),
                 evalue=_fmt_eval(row.get("UNIPROT_evalue")), other=", ".join(other_parts))]


# Order matches the user's own sample-report-before-llm.txt exactly.
_GENERAL_SECTIONS = [rows_rast, rows_pgap, rows_tigrfam, rows_ncbifam, rows_cog, rows_pfam,
                      rows_geneprop, rows_interpro, rows_kegg, rows_eggnog, rows_uniprot]


def rows_tcdb(row: dict) -> list[list[str]]:
    id_ = _clean(row.get("TCDB_subject_id"))
    if not id_:
        return [_row("TCDB")]
    other = _clean(row.get("TCDB_family_description"))
    return [_row("TCDB", id_=id_, desc=_clean(row.get("TCDB_description")),
                 category=_clean(row.get("TCDB_id")), bitscore=_clean(row.get("TCDB_bitscore")),
                 pct_identity=_fmt_pct(row.get("TCDB_percent_identity")),
                 evalue=_fmt_eval(row.get("TCDB_evalue")), other=f"family: {other}" if other else "")]


def rows_merops(row: dict) -> list[list[str]]:
    id_ = _clean(row.get("MEROPS_id"))
    if not id_:
        return [_row("MEROPS")]
    qs, qe = _clean(row.get("MEROPS_query_start")), _clean(row.get("MEROPS_query_end"))
    return [_row("MEROPS", id_=id_, desc=_clean(row.get("MEROPS_description")),
                 domain=f"{qs}-{qe}" if qs and qe else "", category=_clean(row.get("MEROPS_family")),
                 bitscore=_clean(row.get("MEROPS_bitscore")),
                 pct_identity=_fmt_pct(row.get("MEROPS_percent_identity")),
                 evalue=_fmt_eval(row.get("MEROPS_evalue")))]


def rows_dbcan(row: dict) -> list[list[str]]:
    # No Score/Bitscore/%Identity/evalue exist for DBCAN at all (confirmed
    # against the raw table) -- the sample's dash pattern guessed an
    # evalue might exist; it doesn't. EC genuinely IS real here
    # (DBCAN_ec_numbers) even though the sample marked it "-".
    id_ = _clean(row.get("DBCAN_id"))
    if not id_:
        return [_row("DBCAN")]
    n_tools = _clean(row.get("DBCAN_number_of_tools_hit"))
    other = f"{n_tools}/3 detection methods agreed ({_clean(row.get('DBCAN_methods_used'))})" if n_tools else ""
    return [_row("DBCAN", id_=id_, desc=_clean(row.get("DBCAN_description")) or "-",
                 domain=_clean(row.get("DBCAN_hmm_hit")), category=id_,
                 ec=_clean(row.get("DBCAN_ec_numbers")) if _has_value(row.get("DBCAN_ec_numbers")) else "",
                 other=other)]


_SPECIALIZED_SECTIONS = [rows_tcdb, rows_merops, rows_dbcan]
_EVIDENCE_HEADER = ["Tool", "ID", "Description", "Domain", "Category", "EC",
                     "Score", "Bitscore", "%Identity", "E-value", "Other"]


# ---- localization (LLM can read raw topology directly -- not collapsed) ---
# [Tool, Prediction, Domain, Score/Probability, Topology/Cleavage Site, Other]
def rows_signalp6(row: dict) -> list[list[str]]:
    pred = _clean(row.get("SIGNALP6_prediction")) or "OTHER"
    prob_col = {"SP": "SIGNALP6_prob_sp", "LIPO": "SIGNALP6_prob_lipo", "TAT": "SIGNALP6_prob_tat",
                "TATLIPO": "SIGNALP6_prob_tatlipo", "PILIN": "SIGNALP6_prob_pilin"}.get(pred, "SIGNALP6_prob_other")
    other_probs = ",".join(f"{k}:{_clean(row.get(c))}" for k, c in
                            [("sp", "SIGNALP6_prob_sp"), ("lipo", "SIGNALP6_prob_lipo"), ("tat", "SIGNALP6_prob_tat"),
                             ("tatlipo", "SIGNALP6_prob_tatlipo"), ("pilin", "SIGNALP6_prob_pilin"),
                             ("other", "SIGNALP6_prob_other")] if c != prob_col and _clean(row.get(c)))
    return [["SIGNALP6", pred, "", _clean(row.get(prob_col)),
             _clean(row.get("SIGNALP6_cleavage_site")), other_probs]]


def rows_phobius(row: dict) -> list[list[str]]:
    has_sp = _clean(row.get("PHOBIUS_has_signal_peptide")) in ("1", "True", "true")
    n_tm = _clean(row.get("PHOBIUS_n_transmembrane")) or "0"
    pred = "Signal peptide" if has_sp else (f"{n_tm} transmembrane helix(es)" if n_tm not in ("0", "") else "None")
    domain = _clean(row.get("PHOBIUS_segment_label")) or _clean(row.get("PHOBIUS_segment_type"))
    cleave = _clean(row.get("PHOBIUS_signal_cleavage_pos"))
    topo = _clean(row.get("PHOBIUS_topology_short"))
    return [["PHOBIUS", pred, domain, "", f"{topo} (cleavage pos: {cleave})" if cleave else topo, ""]]


def rows_tmbed(row: dict) -> list[list[str]]:
    topo_str = _clean(row.get("TMBED_topology_string"))
    n_seg = _clean(row.get("TMBED_segment_count")) or "1"
    main_topo = _clean(row.get("TMBED_topology"))
    pred = ("Transmembrane" if "transmembrane" in main_topo.lower() else
            "Signal peptide" if "signal" in main_topo.lower() else
            "Non-membrane (cytoplasmic/extracellular)")
    domain = (f"{_clean(row.get('TMBED_segment_start'))}-{_clean(row.get('TMBED_segment_end'))}: {main_topo}"
              if main_topo else "")
    return [["TMBED", pred, domain, "", topo_str, f"{n_seg} segment(s)"]]


def rows_psortb(row: dict) -> list[list[str]]:
    loc = _clean(row.get("PSORTB_localization")) or "Unknown"
    conf = _clean(row.get("PSORTB_is_confident"))
    gram = _clean(row.get("PSORTB_gram_class"))
    return [["PSORTB", loc, "", _clean(row.get("PSORTB_score")), "",
             ", ".join(x for x in [f"confident: {conf}" if conf else "", gram] if x)]]


_LOCALIZATION_SECTIONS = [rows_signalp6, rows_phobius, rows_tmbed, rows_psortb]
_LOCALIZATION_HEADER = ["Tool", "Prediction", "Domain", "Score/Probability", "Topology/Cleavage Site", "Other"]


def evidence_block(section_fns: list, row: dict, header: list[str], indent: str = "") -> str:
    all_rows = [r for fn in section_fns for r in fn(row)]
    return pipe_table(header, all_rows, indent=indent)


# ---- physical neighbor block (peg-number adjacency) -----------------------
def physical_neighbors_block(fid: str, df_indexed: dict, confidence_final: dict,
                              n: int = OPERON_WINDOW) -> str:
    """Compact table of up to n physical genomic neighbors on each side.
    Uses peg-number adjacency as a proxy for chromosomal proximity -- within
    a single replicon, consecutive peg numbers are consecutive on the genome.
    Shown regardless of operon membership so the LLM has genomic context even
    for standalone genes and operon-boundary genes."""
    m = re.match(r'^(fig\|\d+\.\d+\.peg\.)(\d+)$', fid)
    if not m:
        return "Physical neighbor data unavailable (non-peg feature ID)."
    prefix, peg_n = m.group(1), int(m.group(2))
    rows = []
    for offset in list(range(-n, 0)) + list(range(1, n + 1)):
        nbr_fid = f"{prefix}peg.{peg_n + offset}"
        if nbr_fid not in df_indexed:
            continue
        cf_nbr = confidence_final.get(nbr_fid, {})
        label  = (cf_nbr.get("best_consensus_product_descriptor")
                  or cf_nbr.get("best_consensus_product_descriptor")
                  or df_indexed[nbr_fid].get("best_consensus_product_descriptor", "unknown"))
        strand = _clean(df_indexed[nbr_fid].get("RAST_strand", "?"))
        side   = "upstream" if offset < 0 else "downstream"
        rows.append([nbr_fid, f"{offset:+d} ({side})", label[:60], strand])
    if not rows:
        return "No peg-adjacent genes found (contig boundary or single-gene genome)."
    return pipe_table(["Feature ID", "Offset", "Canonical Label", "Strand"], rows)


# ---- operon grouping (sentinel-safe -- see module docstring) --------------
def build_operon_groups(df: pd.DataFrame) -> dict[str, list[str]]:
    cleaned = df["OPERON_id"].astype(str).str.strip()
    has_real = df["OPERON_id"].notna() & (cleaned != "") & ~cleaned.isin(_NO_REAL_OPERON)
    groups: dict[str, list[str]] = {}
    for operon_id, group in df[has_real].groupby("OPERON_id"):
        ordered = group.sort_values("OPERON_gene_position_in_operon")
        groups[operon_id] = ordered["feature_id"].astype(str).tolist()
    return groups


# ---- identity / confidence / fingerprint sections --------------------------
def identity_block(fid: str, row: dict, org_name: str, cf: dict, operon_str: str, indent: str = "") -> str:
    pairs = [
        ("Organism", org_name),
        ("Feature ID", fid),
        ("Gene ID", _clean(row.get("gene_id")) or fid),
        ("Location", f"{_clean(row.get('gene_start'))}..{_clean(row.get('gene_end'))} "
                      f"[{_clean(row.get('RAST_strand'))}]"),
        ("Consensus Product Descriptor", cf.get("best_consensus_product_descriptor", "")),
        ("Product Descriptor Source", cf.get("product_descriptor_source", "")),
        ("Operon", operon_str),
        ("Confidence Score", cf.get("confidence_score", "")),
        ("Confidence Tier", cf.get("confidence_score_tier", "")),
        ("Confidence Flag", cf.get("confidence_flag", "")),
    ]
    return kv(pairs, indent=indent)


def _operon_str(row: dict, full_member_count: int) -> str:
    op_id = _clean(row.get("OPERON_id"))
    if not op_id or op_id in _NO_REAL_OPERON:
        return op_id or "NOT_IN_AN_OPERON"
    pos = _int_str(row.get("OPERON_gene_position_in_operon"))
    return f"{op_id} (position {pos} of {full_member_count})"


def confidence_block(cf: dict) -> str:
    rows = [
        ["C1 Tool Coverage", cf.get("c1_score", ""), cf.get("c1_formula", "")],
        ["C2 Operon Presence", cf.get("c2_score_from_operon_probability", ""), cf.get("c2_formula", "")],
        ["C3 Neighbor Quality", cf.get("c3_score", ""), cf.get("c3_formula", "")],
        ["C4 EC Agreement", cf.get("c4_score", ""), cf.get("c4_formula", "")],
    ]
    lines = [
        "Formula",
        cf.get("confidence_score_formula", "confidence_score = clip([0,1], 0.273 + 0.642*C1 + 0.065*C2 + 0.027*C3 + 0.045*C4)"), "",
        pipe_table(["Component", "Score", "Details"], rows, first_col_w=20), "",
        "Final Score",
        f"Score                   : {cf.get('confidence_score', '')}",
        f"Tier                    : {cf.get('confidence_score_tier', '')}",
        f"Flag                    : {cf.get('confidence_flag', '')}",
    ]
    return "\n".join(lines)


def parse_fingerprint_full(raw: str) -> tuple[str, str, str]:
    h = label = pat = ""
    for part in raw.split(" || "):
        if part.startswith("pattern hash: "):
            h = part[len("pattern hash: "):]
        elif part.startswith("label: "):
            label = part[len("label: "):]
        elif part.startswith("fingerprint: "):
            pat = part[len("fingerprint: "):]
    return h, label, pat


def load_fingerprint_db(db_path: Path, hash_col: str = "fingerprint_hash") -> dict[str, dict]:
    """Load the fingerprint DB once per script run, not per-gene -- an in-memory
    dict avoids re-scanning this ~25k-row file for every gene x fingerprint
    variant (slow on 4000+-gene genomes)."""
    if not db_path.exists():
        return {}
    with open(db_path, newline="") as fh:
        return {row[hash_col]: row for row in csv.DictReader(fh, delimiter="\t") if row.get(hash_col)}


def gene_fingerprint_block(fid: str, fingerprint_full: dict, fingerprint_db: dict[str, dict]) -> str:
    raw = fingerprint_full.get(fid, {}).get("fingerprint", "")
    if not raw:
        return "(no fingerprint -- gene not found in labeled-genes-fingerprint-full.tsv)"
    h, label, pattern = parse_fingerprint_full(raw)
    freq_row = fingerprint_db.get(h, {})
    pairs = [
        ("Hash", h),
        ("Canonical Label", label),
        ("Exact Pattern Frequency", freq_row.get("fingerprint_frequency", "0 (not yet in the cross-genome database)")),
        ("Same Label Frequency", freq_row.get("fingerprint_label_frequency", "0")),
        ("Observed Organisms", freq_row.get("organisms", "(only this organism so far)")),
    ]
    return (kv(pairs) + "\n\nRaw Evidence Pattern\n" + pattern)


_OPERON_FP_VARIANTS = [
    ("Evidence (Ordered)", "operon hash by evidence (ordered): "),
    ("Evidence (Composition)", "operon hash by evidence (composition): "),
    ("Label (Ordered)", "operon hash by label (ordered): "),
    ("Label (Composition)", "operon hash by label (composition): "),
]


def operon_fingerprint_block(fid: str, operon_fp: dict, db_dicts: dict[str, dict]) -> str:
    raw = operon_fp.get(fid, {}).get("operon_fingerprint", "")
    if not raw:
        return "Candidate gene is a standalone (not in an operon) -- no operon fingerprint."
    parts = raw.split(" || ")
    hashes: dict[str, str] = {}
    for variant_label, prefix in _OPERON_FP_VARIANTS:
        for p in parts:
            if p.startswith(prefix):
                hashes[variant_label] = p[len(prefix):]
                break
    rows = []
    for variant_label, _ in _OPERON_FP_VARIANTS:
        h = hashes.get(variant_label, "")
        freq_row = db_dicts[variant_label].get(h, {})
        rows.append([
            variant_label, h,
            freq_row.get("fingerprint_frequency", "0"),
            freq_row.get("fingerprint_label_frequency", "0"),
            freq_row.get("organisms", "(only this organism so far)"),
        ])
    return pipe_table(["Variant", "Hash", "Exact Frequency", "Label Frequency", "Observed Organisms"],
                       rows, first_col_w=24)


# ---- full per-gene document -------------------------------------------------
def operon_member_block(fid: str, df_indexed: dict, confidence_final: dict,
                         org_name: str, operon_groups: dict[str, list[str]]) -> str:
    indent = "    "
    row = df_indexed.get(fid)
    if row is None:
        return f"{indent}[Gene not found in dataset: {fid}]\n"
    cf = confidence_final.get(fid, {})
    op_id = _clean(row.get("OPERON_id"))
    full_member_count = len(operon_groups.get(op_id, [])) or 1
    pos = _int_str(row.get("OPERON_gene_position_in_operon")) or "?"
    operon_str = _operon_str(row, full_member_count)

    lines = [
        f"{indent}OPERON GENE {pos} / {full_member_count}",
        f"{indent}{'-' * 80}", "",
        f"{indent}IDENTITY",
        identity_block(fid, row, org_name, cf, operon_str, indent=indent), "",
        f"{indent}ANNOTATION EVIDENCE", "",
        f"{indent}GENERAL DATABASES",
        evidence_block(_GENERAL_SECTIONS, row, _EVIDENCE_HEADER, indent=indent), "",
        f"{indent}SPECIALIZED DATABASES",
        evidence_block(_SPECIALIZED_SECTIONS, row, _EVIDENCE_HEADER, indent=indent), "",
        f"{indent}LOCALIZATION",
        evidence_block(_LOCALIZATION_SECTIONS, row, _LOCALIZATION_HEADER, indent=indent),
        f"{indent}{'-' * 80}",
    ]
    return "\n".join(lines)


def build_document(fid: str, df_indexed: dict, org_name: str, confidence_final: dict,
                    fingerprint_full: dict, operon_fp: dict, operon_fp_dbs: dict[str, dict],
                    gene_fp_db: dict[str, dict], operon_groups: dict[str, list[str]]) -> str:
    row = df_indexed[fid]
    cf = confidence_final.get(fid, {})

    # Full nesting (every other member's complete evidence block) makes
    # report size scale with operon size -- a 10-gene operon meant every
    # one of its 10 reports repeated ~9x the same content. Capped to the
    # OPERON_WINDOW nearest neighbors on each side instead -- bounds every
    # report to a constant size and matches what actually matters for an
    # operon-coherence check anyway (immediate neighbors), regardless of
    # how large the operon is.
    op_id = _clean(row.get("OPERON_id"))
    other_members, full_member_count = [], 0
    if op_id and op_id not in _NO_REAL_OPERON:
        full_order = operon_groups.get(op_id, [])
        full_member_count = len(full_order)
        if fid in full_order:
            idx = full_order.index(fid)
            other_members = (full_order[max(0, idx - OPERON_WINDOW):idx]
                              + full_order[idx + 1:idx + 1 + OPERON_WINDOW])
        else:
            other_members = [m for m in full_order if m != fid]

    operon_section = "\n".join(
        operon_member_block(m, df_indexed, confidence_final, org_name, operon_groups)
        for m in other_members
    ) if other_members else "Candidate gene is a standalone (not in an operon)."
    operon_truncated = full_member_count - 1 > len(other_members)

    operon_str = _operon_str(row, full_member_count or 1)

    parts = [
        BANNER, "##GENE ANNOTATION REPORT", BANNER, "",
        "IDENTITY",
        identity_block(fid, row, org_name, cf, operon_str), "",
        "##ANNOTATION EVIDENCE", "",
        "#GENERAL DATABASES",
        evidence_block(_GENERAL_SECTIONS, row, _EVIDENCE_HEADER), "",
        "#SPECIALIZED DATABASES",
        evidence_block(_SPECIALIZED_SECTIONS, row, _EVIDENCE_HEADER), "",
        "LOCALIZATION",
        evidence_block(_LOCALIZATION_SECTIONS, row, _LOCALIZATION_HEADER), "",
        "PHYSICAL NEIGHBORS",
        physical_neighbors_block(fid, df_indexed, confidence_final), "",
        "OPERON MEMBERS",
        operon_section, "",
    ]
    if operon_truncated:
        parts += [f"Operon {op_id} has {full_member_count} members total -- showing only the "
                  f"{OPERON_WINDOW} nearest upstream and {OPERON_WINDOW} nearest downstream "
                  f"neighbors above to keep this report's size bounded.", ""]
    parts += [
        "CONFIDENCE SCORE",
        confidence_block(cf), "",
        "GENE FINGERPRINT",
        gene_fingerprint_block(fid, fingerprint_full, gene_fp_db), "",
        "OPERON FINGERPRINT",
        operon_fingerprint_block(fid, operon_fp, operon_fp_dbs), "",
        "SYNTENY",
        "Status",
        "N/A -- synteny is not yet computed for this organism.", "",
        "Comment",
        "Review it the same evidence-grounded way as operon context above once it exists; "
        "until then, do not speculate about it.", "",
        BANNER, "END OF REPORT", BANNER,
    ]
    return "\n".join(parts)


def load_tsv_by_feature(path: Path) -> dict[str, dict]:
    with open(path, newline="") as fh:
        return {r["feature_id"]: r for r in csv.DictReader(fh, delimiter="\t") if r.get("feature_id")}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--consolidated", required=True, help="consolidation/consolidated-merged-all-columns.tsv")
    p.add_argument("--confidence-final", required=True)
    p.add_argument("--fingerprint-full", required=True)
    p.add_argument("--operon-fingerprint", required=True)
    p.add_argument("--fingerprint-database", required=True)
    p.add_argument("--operon-fingerprint-database-evidence-ordered", required=True)
    p.add_argument("--operon-fingerprint-database-evidence-composition", required=True)
    p.add_argument("--operon-fingerprint-database-label-ordered", required=True)
    p.add_argument("--operon-fingerprint-database-label-composition", required=True)
    p.add_argument("--organism-name", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--limit", type=int, default=None)
    return p.parse_args()


def safe_name(s: str) -> str:
    return re.sub(r'[^\w\.\-]', '_', s)


def main() -> None:
    args = parse_args()

    df = pd.read_csv(args.consolidated, sep="\t", low_memory=False)
    df_indexed = {str(r["feature_id"]): r.to_dict() for _, r in df.iterrows()}
    operon_groups = build_operon_groups(df)

    confidence_final = load_tsv_by_feature(Path(args.confidence_final))
    fingerprint_full = load_tsv_by_feature(Path(args.fingerprint_full))
    operon_fp = load_tsv_by_feature(Path(args.operon_fingerprint))

    # Loaded once for the whole run, not once per gene -- see
    # load_fingerprint_db()'s own docstring for why that mattered.
    gene_fp_db = load_fingerprint_db(Path(args.fingerprint_database))
    operon_fp_dbs = {
        "Evidence (Ordered)": load_fingerprint_db(Path(args.operon_fingerprint_database_evidence_ordered)),
        "Evidence (Composition)": load_fingerprint_db(Path(args.operon_fingerprint_database_evidence_composition)),
        "Label (Ordered)": load_fingerprint_db(Path(args.operon_fingerprint_database_label_ordered)),
        "Label (Composition)": load_fingerprint_db(Path(args.operon_fingerprint_database_label_composition)),
    }

    feature_ids = [fid for fid, row in confidence_final.items() if row.get("confidence_score")]
    if args.limit:
        feature_ids = feature_ids[:args.limit]

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    n = 0
    for fid in feature_ids:
        if fid not in df_indexed:
            continue
        doc = build_document(fid, df_indexed, args.organism_name, confidence_final,
                              fingerprint_full, operon_fp, operon_fp_dbs,
                              gene_fp_db, operon_groups)
        (output_dir / f"{safe_name(fid)}.txt").write_text(doc)
        n += 1

    print(f"[build-gene-report] prepared {n} gene reports -> {output_dir}")


if __name__ == "__main__":
    main()
