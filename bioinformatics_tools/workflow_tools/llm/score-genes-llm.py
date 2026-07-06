#!/usr/bin/env python3
"""score-genes-llm.py — margie_sb phase15 (llm), step 3 (GPU): enhanced LLM
confidence scoring layer on top of the mechanical phase11 scores.

The LLM receives the full gene annotation report (evidence, localization,
operon context, mechanical C1-C4, fingerprints) and does four things the
mechanical pipeline cannot:

  C3 re-assessment  — semantic pathway coherence: do operon neighbors
                      actually corroborate this gene's function?
  C4 re-assessment  — EC agreement with granularity awareness: a conflict
                      between 2.7.1.2 and 2.7.1.- is NOT a real conflict;
                      only contradictory reaction types/substrates are.
  C5 scoring        — conflict/ambiguity: pairwise reasoning across all
                      database hits to detect genuine functional conflicts
                      vs. nomenclature/granularity differences.
  Topology check    — does the canonical label's expected
                      compartment/topology match SignalP/Phobius/TMbed/PSORTb?

C1 (tool coverage fraction) and C2 (operon reliability geometric mean) are
accepted from the mechanical pipeline unchanged — they are pure arithmetic.

Final LLM confidence:
  formula_llm = clip([0,1],
      0.111 + 0.741*mech_C1 + 0.063*mech_C2
            + 0.111*llm_C3  + 0.068*llm_C4  - 0.056*llm_C5
            + topology_adj)
  topology_adj = -0.05 if TOPOLOGY is "inconsistent", else 0.

The LLM also outputs a holistic LLM_CONFIDENCE directly (0-1). Both the
formula-based and holistic scores are written to the summary TSV; the
holistic score is the one used downstream as the final confidence.

Usage:
  python score-genes-llm.py \\
    --trained-model  <path to model dir, base or LoRA adapter> \\
    --prepared-dir   <phase14 build-gene-report.py output dir> \\
    --reports-dir    <dir for final per-gene .txt reports> \\
    --summary        <output summary TSV path> \\
    [--max-tokens 100000] [--limit N]
"""
from __future__ import annotations

import argparse
import csv
import json
import logging
import re
import sys
from pathlib import Path

LOGGER = logging.getLogger("score_genes_llm")

H = "=" * 100

# ---------------------------------------------------------------------------
# Calibration coefficients (manual formula, from scoring-development study)
# C1 and C2: mechanical values passed through unchanged
# C3, C4, C5: LLM-assessed values replace the mechanical ones
# ---------------------------------------------------------------------------
_INTERCEPT = 0.111
_W = {"c1": 0.741, "c2": 0.063, "c3": 0.111, "c4": 0.068, "c5": -0.056}
_TOPOLOGY_PENALTY = -0.05   # applied when TOPOLOGY == "inconsistent"

SYSTEM_PROMPT = """\
You are an expert microbial genome annotator reviewing a single gene's \
annotation evidence report produced by the MARGIE(SB) pipeline.

The report contains:
  - The gene's canonical label and its mechanical confidence score (C1-C4)
  - Full database evidence from up to 15 annotation tools
  - Localization predictions (SignalP6, Phobius, TMBED, PSORTb)
  - Physical genomic neighbors (±3 peg-adjacent genes: label, strand, offset)
  - Operon members (if the gene is predicted to be in an operon)
  - Gene and operon fingerprints

The pipeline has already computed mechanical scores for C1 and C2.
  C1 (Tool Coverage) — fraction of tools that returned an informative hit.
      Accept this value as-is. It is simple arithmetic.
  C2 (Operon Reliability) — geometric mean of operon-member probabilities.
      Accept this value as-is. It is simple arithmetic.

Your job is to assess the three remaining components that the mechanical
pipeline cannot evaluate accurately:

  C3 (Neighborhood Coherence — two structured questions, NOT a raw float):
      Look at the PHYSICAL NEIGHBORS section (±3 peg-adjacent genes,
      regardless of operon boundaries). Answer two questions independently:

      Q1 — NEIGHBORHOOD_FIT: Does this gene's canonical label belong to the
           same functional theme as any of the physical neighbors?
           strong   — label directly identifies the gene as a member of a
                      recognizable functional unit present among the neighbors
                      (e.g., neighbors are all ABC transporter components and
                      this gene's label is also an ABC transporter component;
                      or neighbors form a known biosynthesis pathway and this
                      gene's label names a role in that pathway)
           moderate — label is functionally consistent with the dominant
                      neighborhood theme but does not identify a specific
                      shared functional unit
           none     — label shares no functional theme with neighbors, OR the
                      gene has fewer than 2 identifiable neighbors, OR all
                      neighbors are hypothetical proteins (absence of fit is
                      neutral, not negative)

      Q2 — NEIGHBORHOOD_CONTRADICTION: Does this gene's label functionally
           contradict its neighbors?
           no  — no contradiction (the default; "different function" is NOT
                 contradiction; hypothetical neighbors are always "no")
           yes — label implies a function that is clearly incongruous with
                 the neighborhood's dominant functional theme (e.g., a gene
                 labeled "transcriptional repressor" whose 6 neighbors are
                 all ribosomal proteins — not just different, but jarring)

      C3 is computed mechanically from your two answers — do NOT output a
      raw LLM_C3 float. Output NEIGHBORHOOD_FIT and NEIGHBORHOOD_CONTRADICTION
      in the SCORES block instead. The formula is:
        fit_bonus     = strong:+0.4, moderate:+0.2, none:0.0
        contradiction = yes:−0.2, no:0.0
        C3 = clip(0.5 + fit_bonus + contradiction, min=0.3, max=0.9)

  C4 (EC Agreement, score 0.00–1.00):
      CRITICAL: Only examine EC numbers reported for the CANDIDATE GENE in the
      ##ANNOTATION EVIDENCE section (#GENERAL DATABASES + #SPECIALIZED DATABASES).
      Do NOT use EC numbers from the OPERON MEMBERS sections — those belong to
      neighbouring genes and must NOT influence C4 for the candidate gene.
      Among the candidate gene's own EC numbers:
        Different levels of EC specificity are NOT conflicts.
          2.7.1.2  (glucokinase) vs  2.7.1.-  (sugar kinase) → same enzyme,
          different resolution → NOT a conflict, treat as near-agreement.
          2.7.1.2  vs  1.1.1.27 → different reaction class → genuine conflict.
      Score 1.0 for full consensus. Score near 0 only for genuine reaction-
      class conflicts. Score 0.5 when the candidate gene has no EC assignment.

  C5 (Conflict/Ambiguity, score 0.00–1.00):
      CRITICAL: Only consider database hits for the CANDIDATE GENE (from the
      ##ANNOTATION EVIDENCE section — #GENERAL DATABASES + #SPECIALIZED DATABASES).
      Do NOT include any hits from the OPERON MEMBERS sections.
      Read the candidate gene's informative hits as a set. For every pair, decide:
        Agree    — same function (different names/IDs are fine)
        Ambiguous — one is broader/more general than the other
        Conflict  — genuinely different functions that cannot both be true
      C5 = (n_conflict + n_ambiguous) / n_pairs_compared
      Score 0.0 = perfect consensus. Score 1.0 = maximum disagreement.
      If fewer than two informative hits exist, score 0.0.

  Topology & Localization check:
      Compare the canonical label's implied cellular location/topology with
      the localization predictions (SignalP6, Phobius, TMbed, PSORTb).
      Examples: a "secreted protease" should have a signal peptide; an
      "inner membrane transporter" should have transmembrane helices; a
      "cytoplasmic enzyme" should have no predicted signal/TM topology.
      Classify as: consistent | inconsistent | not_applicable
      not_applicable: label gives no localization expectation, or all tools
      returned "no prediction".

  Specialized Database Cross-Check:
      Examine hits from TIGRFAM, HAMAP (often under UNIPROT with HPA/MF
      accessions), CDD, and PSORTb. These are curated, high-specificity
      databases with family-level precision.
      If any returned a functional hit, does its description agree with the
      canonical label?
        agrees            — all specialized DB hits are consistent with the label
        disagrees         — at least one specialized DB hit contradicts the label
        partial           — some agree, some disagree
        no_specialized_hit — none of TIGRFAM/HAMAP/CDD returned a hit

  Operon Coherence:
      Based on your C3 pathway coherence assessment, classify this gene's
      operon context:
        yes            — operon neighbors functionally corroborate the label
        no             — operon neighbors contradict or are inconsistent with label
        not_applicable — standalone gene, no operon, or all neighbors hypothetical

Write your analysis first, using exactly these section headers:

## DATABASE COVERAGE
[2-3 sentences: which databases returned informative hits, quality and \
specificity of hits, how well C1 reflects actual functional certainty]

## CONFLICTS AND AMBIGUITY
[2-3 sentences: read the CANDIDATE GENE'S OWN database hits together \
(##ANNOTATION EVIDENCE only — not operon neighbours) — do they agree? \
any genuine functional contradictions vs naming/granularity differences? derive C5]

## NEIGHBORHOOD COHERENCE
[2-3 sentences: look at the PHYSICAL NEIGHBORS section. What is the dominant \
functional theme of the ±3 genomic neighbors? State your NEIGHBORHOOD_FIT \
answer (strong/moderate/none) and your NEIGHBORHOOD_CONTRADICTION answer \
(no/yes) explicitly, with one sentence of reasoning for each.]

## EC ASSESSMENT
[2-3 sentences: which tools reported EC numbers for the CANDIDATE GENE ONLY \
(##ANNOTATION EVIDENCE section — not operon neighbours)? do they agree? are \
any apparent conflicts just different levels of EC specificity? derive C4]

## TOPOLOGY AND LOCALIZATION
[1-2 sentences: does the canonical label's expected cellular location match \
SignalP/Phobius/TMbed/PSORTb predictions? note any inconsistency explicitly]

Then end with EXACTLY this block (no extra text after it):

## SCORES
NEIGHBORHOOD_FIT: <strong | moderate | none>
NEIGHBORHOOD_CONTRADICTION: <no | yes>
LLM_C4: <0.00–1.00>
LLM_C5: <0.00–1.00>
TOPOLOGY: <consistent | inconsistent | not_applicable>
SPECIALIZED_DB_AGREEMENT: <agrees | disagrees | partial | no_specialized_hit>
OPERON_COHERENT: <yes | no | not_applicable>
LLM_CONFIDENCE: <0.00–1.00>
SUMMARY: <one sentence, plain language, suitable for a reviewer scanning \
hundreds of genes>

Rules for LLM_CONFIDENCE:
  - If the mechanical score already reflects the evidence accurately, copy it.
  - If your C3/C4/C5 assessment reveals that the mechanical score is too \
high or too low, adjust it accordingly and state why in the relevant section.
  - This is your holistic judgment of annotation reliability, not a formula.
"""

BOS   = "<|begin_of_text|>"
SYS_S = "<|start_header_id|>system<|end_header_id|>\n\n"
SYS_E = "<|eot_id|>"
USR_S = "<|start_header_id|>user<|end_header_id|>\n\n"
USR_E = "<|eot_id|>"
ASS_S = "<|start_header_id|>assistant<|end_header_id|>\n\n"


def build_prompt(review_document: str) -> str:
    return (
        BOS
        + SYS_S + SYSTEM_PROMPT + SYS_E
        + USR_S + review_document + USR_E
        + ASS_S
    )


def _detect_backend() -> str:
    try:
        import mlx_lm  # noqa: F401
        return "mlx"
    except ImportError:
        pass
    try:
        import torch  # noqa: F401
        import transformers  # noqa: F401
        return "hf"
    except ImportError:
        pass
    return "none"


_BACKEND = _detect_backend()


def load_model(model_path: Path):
    if not (model_path / "config.json").exists() and not (model_path / "adapter_config.json").exists():
        LOGGER.error(
            f"No config.json or adapter_config.json under {model_path} -- "
            "no model weights found."
        )
        sys.exit(1)

    LOGGER.info(f"Loading model from {model_path} ...")
    if _BACKEND == "mlx":
        from mlx_lm import load
        model, tokenizer = load(str(model_path))
        LOGGER.info("Model loaded via mlx_lm (Apple Silicon).")
        return model, tokenizer
    elif _BACKEND == "hf":
        import torch
        from transformers import AutoModelForCausalLM, AutoTokenizer
        adapter_cfg = model_path / "adapter_config.json"
        if adapter_cfg.exists():
            from peft import PeftModel
            base_name = json.loads(adapter_cfg.read_text()).get("base_model_name_or_path")
            LOGGER.info(f"  Detected LoRA adapter. Base model: {base_name}")
            tokenizer = AutoTokenizer.from_pretrained(base_name)
            tokenizer.pad_token = tokenizer.eos_token
            base = AutoModelForCausalLM.from_pretrained(
                base_name, torch_dtype=torch.bfloat16, device_map="auto",
            )
            model = PeftModel.from_pretrained(base, str(model_path))
            model.eval()
        else:
            tokenizer = AutoTokenizer.from_pretrained(str(model_path))
            tokenizer.pad_token = tokenizer.eos_token
            model = AutoModelForCausalLM.from_pretrained(
                str(model_path), torch_dtype=torch.bfloat16, device_map="auto",
            )
            model.eval()
        LOGGER.info(f"Model loaded via transformers on {next(model.parameters()).device}.")
        return model, tokenizer
    else:
        LOGGER.error(
            "No inference backend found. Install mlx-lm (macOS) or "
            "torch+transformers (Linux/HPC)."
        )
        sys.exit(1)


def generate_response(model, tokenizer, prompt: str, max_tokens: int) -> str:
    if _BACKEND == "mlx":
        from mlx_lm import generate
        output = generate(model, tokenizer, prompt=prompt, max_tokens=max_tokens, verbose=False)
        if output.startswith(prompt):
            output = output[len(prompt):]
        text = output.strip()
    else:
        import torch
        inputs = tokenizer(prompt, return_tensors="pt").to(model.device)
        with torch.no_grad():
            output_ids = model.generate(
                **inputs, max_new_tokens=max_tokens, do_sample=False,
                temperature=None, top_p=None, pad_token_id=tokenizer.eos_token_id,
            )
        text = tokenizer.decode(
            output_ids[0][inputs["input_ids"].shape[1]:], skip_special_tokens=True,
        ).strip()
    return re.sub(r'\s*\bSTOP\b\s*$', '', text, flags=re.IGNORECASE).strip()


def generate_batch(model, tokenizer, prompts: list[str], max_tokens: int) -> list[str]:
    """Generate responses for a batch of prompts simultaneously on GPU (HF only).
    Falls back to sequential on MLX."""
    if _BACKEND != "hf" or len(prompts) == 1:
        return [generate_response(model, tokenizer, p, max_tokens) for p in prompts]

    import torch
    # Left-pad so all generated tokens land after the (right-side) input end
    orig_padding_side = tokenizer.padding_side
    tokenizer.padding_side = "left"
    try:
        inputs = tokenizer(
            prompts, return_tensors="pt", padding=True, truncation=False,
        ).to(model.device)
        input_len = inputs["input_ids"].shape[1]
        with torch.no_grad():
            output_ids = model.generate(
                **inputs, max_new_tokens=max_tokens, do_sample=False,
                temperature=None, top_p=None, pad_token_id=tokenizer.eos_token_id,
            )
        texts = tokenizer.batch_decode(
            output_ids[:, input_len:], skip_special_tokens=True,
        )
    finally:
        tokenizer.padding_side = orig_padding_side

    return [re.sub(r'\s*\bSTOP\b\s*$', '', t, flags=re.IGNORECASE).strip() for t in texts]


def _float_field(text: str, key: str) -> float | None:
    m = re.search(rf'^{re.escape(key)}\s*:\s*([0-9]+(?:\.[0-9]+)?)', text,
                  re.IGNORECASE | re.MULTILINE)
    if not m:
        return None
    v = float(m.group(1))
    return max(0.0, min(1.0, v))


def _str_field(text: str, key: str) -> str:
    m = re.search(rf'^{re.escape(key)}\s*:\s*(.+?)(?:\n|\Z)', text,
                  re.IGNORECASE | re.MULTILINE)
    return m.group(1).strip() if m else ""


def _section_text(text: str, header: str) -> str:
    """Extract the body of a ## HEADER section."""
    m = re.search(
        rf'##\s*{re.escape(header)}\s*\n(.*?)(?=\n##\s|\Z)',
        text, re.IGNORECASE | re.DOTALL,
    )
    return m.group(1).strip() if m else ""


def _compute_verdict(llm_confidence: float | None, mech_score: float | None) -> str:
    if llm_confidence is None or mech_score is None:
        return "UNPARSEABLE"
    delta = abs(llm_confidence - mech_score)
    if delta <= 0.10:
        return "AGREES"
    if delta <= 0.25:
        return "PARTIAL"
    return "DISAGREES"


def parse_llm_output(response: str, mech_score: float | None = None) -> dict:
    scores_block = _section_text(response, "SCORES")
    if not scores_block:
        scores_block = response

    nbr_fit     = _str_field(scores_block, "NEIGHBORHOOD_FIT").lower()
    nbr_cont    = _str_field(scores_block, "NEIGHBORHOOD_CONTRADICTION").lower()
    llm_c4      = _float_field(scores_block, "LLM_C4")
    llm_c5      = _float_field(scores_block, "LLM_C5")
    topology    = _str_field(scores_block, "TOPOLOGY").lower()
    spec_agree  = _str_field(scores_block, "SPECIALIZED_DB_AGREEMENT").lower()
    operon_coh  = _str_field(scores_block, "OPERON_COHERENT").lower()
    llm_conf    = _float_field(scores_block, "LLM_CONFIDENCE")
    summary     = _str_field(scores_block, "SUMMARY")

    # Compute C3 mechanically from the two categorical answers
    _FIT_BONUS   = {"strong": 0.4, "moderate": 0.2, "none": 0.0}
    _CONT_PENALTY = {"yes": -0.2, "no": 0.0}
    fit_bonus     = _FIT_BONUS.get(nbr_fit, 0.0)
    cont_penalty  = _CONT_PENALTY.get(nbr_cont, 0.0)
    llm_c3 = max(0.3, min(0.9, 0.5 + fit_bonus + cont_penalty))

    _valid_spec = {"agrees", "disagrees", "partial", "no_specialized_hit"}
    _valid_coh  = {"yes", "no", "not_applicable"}

    return {
        "neighborhood_fit":          nbr_fit if nbr_fit in ("strong","moderate","none") else "unparseable",
        "neighborhood_contradiction": nbr_cont if nbr_cont in ("yes","no") else "unparseable",
        "llm_c3":       llm_c3,
        "llm_c4":       llm_c4,
        "llm_c5":       llm_c5,
        "topology":     topology if topology in ("consistent","inconsistent","not_applicable") else "unparseable",
        "specialized_db_agreement": spec_agree if spec_agree in _valid_spec else "unparseable",
        "operon_coherent":          operon_coh if operon_coh in _valid_coh  else "unparseable",
        "llm_confidence": llm_conf,
        "verdict":      _compute_verdict(llm_conf, mech_score),
        "summary":      summary,
        "db_coverage":      _section_text(response, "DATABASE COVERAGE"),
        "conflicts":        _section_text(response, "CONFLICTS AND AMBIGUITY"),
        "neighborhood_ev":  _section_text(response, "NEIGHBORHOOD COHERENCE"),
        "ec_assessment":    _section_text(response, "EC ASSESSMENT"),
        "topology_text":    _section_text(response, "TOPOLOGY AND LOCALIZATION"),
    }


def compute_formula_score(mech_c1: float, mech_c2: float,
                           llm_c3: float | None, llm_c4: float | None,
                           llm_c5: float | None, topology: str) -> float | None:
    """Manual-calibration formula with LLM-assessed C3/C4/C5."""
    if any(v is None for v in (llm_c3, llm_c4, llm_c5)):
        return None
    adj = _TOPOLOGY_PENALTY if topology == "inconsistent" else 0.0
    raw = (
        _INTERCEPT
        + _W["c1"] * mech_c1
        + _W["c2"] * mech_c2
        + _W["c3"] * llm_c3
        + _W["c4"] * llm_c4
        + _W["c5"] * llm_c5
        + adj
    )
    return max(0.0, min(1.0, raw))


def _extract_field(doc: str, label: str) -> str:
    m = re.search(rf'{re.escape(label)}\s*:\s*(.+?)(?:\n|\Z)', doc)
    return m.group(1).strip() if m else ""


def _safe_float_field(doc: str, label: str) -> float:
    try:
        return float(_extract_field(doc, label))
    except ValueError:
        return 0.5


# ─────────────────────────────────────────────────────────────────────────────
# Report field parser — extracts every traceability column from the prepared
# gene annotation report produced by build-gene-report.py (phase 14).
# ─────────────────────────────────────────────────────────────────────────────

_GENERAL_TOOLS = ["RAST", "PGAP", "TIGRFAM", "NCBIFAM", "COG", "PFAM",
                   "GENEPROP", "INTERPRO", "KEGG", "EGGNOG", "UNIPROT"]
_SPECIAL_TOOLS = ["TCDB", "MEROPS", "DBCAN"]
_LOCAL_TOOLS   = ["SIGNALP6", "PHOBIUS", "TMBED", "PSORTB"]


def _table_rows(block: str) -> dict:
    """Parse all pipe-delimited rows in a text block.
    Returns {TOOL: [[col1, col2, ...], ...]} (tool column excluded)."""
    rows: dict = {}
    for line in block.split("\n"):
        if "|" not in line:
            continue
        cols = [c.strip() for c in line.split("|")]
        if len(cols) < 2:
            continue
        tool = cols[0]
        if not tool or tool == "Tool" or re.match(r"^[-\s]+$", tool):
            continue
        rows.setdefault(tool, []).append(cols[1:])
    return rows


def _join_vals(vals: list) -> str:
    seen: list = []
    for v in vals:
        if v and v != "-" and v not in seen:
            seen.append(v)
    return "; ".join(seen)


def _component_score(block: str, label: str) -> str:
    m = re.search(rf"^{re.escape(label)}\s*\|\s*([0-9.]+)", block, re.MULTILINE)
    return m.group(1) if m else ""


def parse_report_fields(doc: str) -> dict:
    """Extract every traceability field from a prepared gene annotation report."""
    f: dict = {}

    # ── IDENTITY ──────────────────────────────────────────────────────────────
    f["organism"]                   = _extract_field(doc, "Organism")
    f["domain"]                     = ""   # comes from GTDBTk, not in report
    f["feature_id"]                 = _extract_field(doc, "Feature ID")
    f["gene_id"]                    = _extract_field(doc, "Gene ID")
    f["canonical_concordant_label"] = _extract_field(doc, "Canonical Label")
    f["canonical_source"]           = _extract_field(doc, "Canonical Source")

    loc = _extract_field(doc, "Location")
    lm  = re.match(r"(\d+)\.\.(\d+)\s*\[([+\-])\]", loc)
    f["start"]  = lm.group(1) if lm else ""
    f["end"]    = lm.group(2) if lm else ""
    f["strand"] = lm.group(3) if lm else ""

    # ── Only parse the main-gene section (before OPERON MEMBERS) ─────────────
    om_m = re.search(r"^OPERON MEMBERS", doc, re.MULTILINE)
    pre  = doc[:om_m.start()] if om_m else doc

    gdb_m  = re.search(r"#GENERAL DATABASES",    pre)
    sdb_m  = re.search(r"#SPECIALIZED DATABASES", pre)
    loc_m  = re.search(r"^LOCALIZATION",           pre, re.MULTILINE)
    cs_m   = re.search(r"^CONFIDENCE SCORE",       pre, re.MULTILINE)

    gdb_block = pre[gdb_m.end()  : sdb_m.start()] if gdb_m and sdb_m else ""
    sdb_block = pre[sdb_m.end()  : loc_m.start()]  if sdb_m and loc_m else ""
    loc_block = pre[loc_m.end()  : cs_m.start()]   if loc_m and cs_m  else (pre[loc_m.end():] if loc_m else "")

    gdb_rows = _table_rows(gdb_block)
    sdb_rows = _table_rows(sdb_block)
    loc_rows = _table_rows(loc_block)

    # ── General annotation databases ──────────────────────────────────────────
    for tool in _GENERAL_TOOLS:
        key  = tool.lower()
        rows = gdb_rows.get(tool, [])
        ids  = _join_vals([r[0] for r in rows])
        desc = _join_vals([r[1] for r in rows if len(r) > 1])
        if tool == "RAST":
            f["rast_description"] = desc
        else:
            f[f"{key}_id"]          = ids
            f[f"{key}_description"] = desc

    # ── Specialized databases ──────────────────────────────────────────────────
    for tool in _SPECIAL_TOOLS:
        key  = tool.lower()
        rows = sdb_rows.get(tool, [])
        f[f"{key}_id"]          = _join_vals([r[0] for r in rows])
        f[f"{key}_description"] = _join_vals([r[1] for r in rows if len(r) > 1])

    # ── Localization ──────────────────────────────────────────────────────────
    # Table cols after tool: Prediction | Domain | Score/Prob | Topology/CS | Other
    def _loc_col(tool: str, idx: int) -> str:
        rows = loc_rows.get(tool, [])
        v = rows[0][idx] if rows and len(rows[0]) > idx else ""
        return "" if v == "-" else v

    f["signalp6_prediction"] = _loc_col("SIGNALP6", 0)
    f["signalp6_score"]      = _loc_col("SIGNALP6", 2)
    f["phobius_prediction"]  = _loc_col("PHOBIUS",  0)
    f["phobius_topology"]    = _loc_col("PHOBIUS",  3)
    f["tmbed_prediction"]    = _loc_col("TMBED",    0)
    f["psortb_prediction"]   = _loc_col("PSORTB",   0)
    f["psortb_score"]        = _loc_col("PSORTB",   2)

    # ── CONFIDENCE SCORE section (after OPERON MEMBERS in the full doc) ───────
    cs_full = re.search(r"^CONFIDENCE SCORE", doc, re.MULTILINE)
    cs_block = doc[cs_full.start():] if cs_full else ""

    f["mechanical_c1"]   = _component_score(cs_block, "C1 Tool Coverage")
    f["mechanical_c2"]   = _component_score(cs_block, "C2 Operon Presence")
    f["mechanical_c3"]   = _component_score(cs_block, "C3 Neighbor Quality")
    f["mechanical_c4"]   = _component_score(cs_block, "C4 EC Agreement")
    f["mechanical_score"] = _extract_field(cs_block, "Score")
    f["mechanical_tier"]  = _extract_field(cs_block, "Tier")
    f["mechanical_flag"]  = _extract_field(cs_block, "Flag")

    # ── GENE FINGERPRINT section ──────────────────────────────────────────────
    fp_m   = re.search(r"^GENE FINGERPRINT",   doc, re.MULTILINE)
    opfp_m = re.search(r"^OPERON FINGERPRINT", doc, re.MULTILINE)
    fp_block = doc[fp_m.start(): opfp_m.start()] if fp_m and opfp_m else ""

    f["fingerprint_hash"]                       = _extract_field(fp_block, "Hash")
    f["fingerprint_hash_frequency_in_database"] = _extract_field(fp_block, "Exact Pattern Frequency")
    f["label_frequency_in_database"]            = _extract_field(fp_block, "Same Label Frequency")

    raw_m = re.search(r"Raw Evidence Pattern\s*\n(.+?)(?:\n\n|\Z)", fp_block, re.DOTALL)
    f["associated_fingerprint"]                       = raw_m.group(1).strip() if raw_m else ""
    f["associated_fingerprint_frequency_in_database"] = f["fingerprint_hash_frequency_in_database"]

    return f


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--trained-model", required=True)
    p.add_argument("--prepared-dir",  required=True)
    p.add_argument("--reports-dir",   required=True)
    p.add_argument("--summary",       required=True)
    p.add_argument("--max-tokens",  type=int, default=800)
    p.add_argument("--limit",       type=int, default=None)
    p.add_argument("--batch-size",  type=int, default=1,
                   help="Number of gene reports to score in a single GPU call (default: 1)")
    return p.parse_args()


SUMMARY_COLUMNS = [
    # Identity
    "organism", "domain", "feature_id", "gene_id", "start", "end", "strand",
    # Label
    "canonical_concordant_label", "canonical_source",
    # Fingerprint
    "label_frequency_in_database",
    "fingerprint_hash", "fingerprint_hash_frequency_in_database",
    "associated_fingerprint", "associated_fingerprint_frequency_in_database",
    # General annotation databases
    "rast_description",
    "pgap_id", "pgap_description",
    "tigrfam_id", "tigrfam_description",
    "ncbifam_id", "ncbifam_description",
    "cog_id", "cog_description",
    "pfam_id", "pfam_description",
    "geneprop_id", "geneprop_description",
    "interpro_id", "interpro_description",
    "kegg_id", "kegg_description",
    "eggnog_id", "eggnog_description",
    "uniprot_id", "uniprot_description",
    # Specialized databases
    "tcdb_id", "tcdb_description",
    "merops_id", "merops_description",
    "dbcan_id", "dbcan_description",
    # Localization
    "signalp6_prediction", "signalp6_score",
    "phobius_prediction", "phobius_topology",
    "tmbed_prediction",
    "psortb_prediction", "psortb_score",
    # Mechanical assessment
    "mechanical_c1", "mechanical_c2", "mechanical_c3", "mechanical_c4",
    "mechanical_score", "mechanical_tier", "mechanical_flag",
    # LLM assessment
    "neighborhood_fit", "neighborhood_contradiction",
    "llm_c3", "llm_c4", "llm_c5", "topology",
    "specialized_db_agreement",
    "operon_coherent",
    "formula_llm_score",
    "llm_assessment_score",
    "llm_agreement_with_label",
    "llm_confidence_on_label",
    "llm_reasoning_text",
]


def main() -> None:
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s  %(levelname)s  %(message)s")
    args = parse_args()

    prepared_dir   = Path(args.prepared_dir)
    prepared_files = sorted(prepared_dir.glob("*.txt"))
    if args.limit:
        prepared_files = prepared_files[:args.limit]
    LOGGER.info(f"{len(prepared_files)} prepared documents in {prepared_dir}")

    model, tokenizer = load_model(Path(args.trained_model))

    reports_dir = Path(args.reports_dir)
    reports_dir.mkdir(parents=True, exist_ok=True)
    summary_path = Path(args.summary)
    summary_path.parent.mkdir(parents=True, exist_ok=True)

    verdict_counts: dict[str, int] = {}
    n_done = 0
    batch_size = max(1, args.batch_size)
    LOGGER.info(f"Batch size: {batch_size} ({'HF batched' if _BACKEND == 'hf' and batch_size > 1 else 'sequential'})")

    with open(summary_path, "w", newline="") as out_fh:
        writer = csv.DictWriter(out_fh, fieldnames=SUMMARY_COLUMNS, delimiter="\t")
        writer.writeheader()

        for batch_start in range(0, len(prepared_files), batch_size):
            batch_paths = prepared_files[batch_start : batch_start + batch_size]
            review_docs = [p.read_text() for p in batch_paths]
            prompts     = [build_prompt(doc) for doc in review_docs]

            LOGGER.info(f"Scoring batch {batch_start // batch_size + 1}: "
                        f"genes {batch_start + 1}–{batch_start + len(batch_paths)} "
                        f"of {len(prepared_files)}")
            llm_texts = generate_batch(model, tokenizer, prompts, args.max_tokens)

            for prepared_path, review_doc, llm_text in zip(batch_paths, review_docs, llm_texts):
                mech_score_str = _extract_field(review_doc, "Confidence Score")
                mech_tier      = _extract_field(review_doc, "Confidence Tier")
                try:
                    mech_score_f = float(mech_score_str)
                except (ValueError, TypeError):
                    mech_score_f = None

                # Pull mechanical C1/C2 from the report for the formula calculation
                mech_c1 = _safe_float_field(review_doc, "C1 Tool Coverage")
                mech_c2 = _safe_float_field(review_doc, "C2 Operon Presence")

                parsed = parse_llm_output(llm_text, mech_score=mech_score_f)

                formula_score = compute_formula_score(
                    mech_c1, mech_c2,
                    parsed["llm_c3"], parsed["llm_c4"], parsed["llm_c5"],
                    parsed["topology"],
                )

                mech_score = mech_score_str

                # Build the LLM section appended to the report file
                llm_section = (
                    f"\n{H}\n"
                    f"LLM CONFIDENCE ASSESSMENT\n"
                    f"{H}\n\n"
                    f"## DATABASE COVERAGE\n{parsed['db_coverage']}\n\n"
                    f"## CONFLICTS AND AMBIGUITY\n{parsed['conflicts']}\n\n"
                    f"## NEIGHBORHOOD COHERENCE\n{parsed['neighborhood_ev']}\n\n"
                    f"## EC ASSESSMENT\n{parsed['ec_assessment']}\n\n"
                    f"## TOPOLOGY AND LOCALIZATION\n{parsed['topology_text']}\n\n"
                    f"## SCORES\n"
                    f"NEIGHBORHOOD_FIT         : {parsed['neighborhood_fit']}\n"
                    f"NEIGHBORHOOD_CONTRADICTION: {parsed['neighborhood_contradiction']}\n"
                    "LLM_C3 (computed)        : " + (f"{parsed['llm_c3']:.4f}" if parsed['llm_c3'] is not None else 'unparseable') + "\n"
                    "LLM_C4                   : " + (f"{parsed['llm_c4']}" if parsed['llm_c4'] is not None else 'unparseable') + "\n"
                    "LLM_C5                   : " + (f"{parsed['llm_c5']}" if parsed['llm_c5'] is not None else 'unparseable') + "\n"
                    f"TOPOLOGY                 : {parsed['topology']}\n"
                    f"SPECIALIZED_DB_AGREEMENT : {parsed['specialized_db_agreement']}\n"
                    f"OPERON_COHERENT          : {parsed['operon_coherent']}\n"
                    "Formula LLM Score        : " + (f"{formula_score:.4f}" if formula_score is not None else 'unparseable') + "\n"
                    "LLM_CONFIDENCE           : " + (f"{parsed['llm_confidence']:.4f}" if parsed['llm_confidence'] is not None else 'unparseable') + "\n"
                    f"VERDICT                  : {parsed['verdict']}\n"
                    f"SUMMARY                  : {parsed['summary']}\n"
                )

                full_report = review_doc + llm_section
                (reports_dir / prepared_path.name).write_text(full_report)

                rf = parse_report_fields(review_doc)

                llm_reasoning = " | ".join(filter(None, [
                    parsed["db_coverage"].replace("\n", " ").strip(),
                    parsed["conflicts"].replace("\n", " ").strip(),
                    parsed["neighborhood_ev"].replace("\n", " ").strip(),
                    parsed["ec_assessment"].replace("\n", " ").strip(),
                    parsed["topology_text"].replace("\n", " ").strip(),
                    parsed["summary"].strip() if parsed["summary"] else "",
                ]))

                writer.writerow({
                    # Identity
                    "organism":    rf["organism"],
                    "domain":      rf["domain"],
                    "feature_id":  rf["feature_id"],
                    "gene_id":     rf["gene_id"],
                    "start":       rf["start"],
                    "end":         rf["end"],
                    "strand":      rf["strand"],
                    # Label
                    "canonical_concordant_label": rf["canonical_concordant_label"],
                    "canonical_source":           rf["canonical_source"],
                    # Fingerprint
                    "label_frequency_in_database":                rf["label_frequency_in_database"],
                    "fingerprint_hash":                           rf["fingerprint_hash"],
                    "fingerprint_hash_frequency_in_database":     rf["fingerprint_hash_frequency_in_database"],
                    "associated_fingerprint":                     rf["associated_fingerprint"],
                    "associated_fingerprint_frequency_in_database": rf["associated_fingerprint_frequency_in_database"],
                    # General annotation databases
                    "rast_description":    rf.get("rast_description", ""),
                    "pgap_id":             rf.get("pgap_id", ""),
                    "pgap_description":    rf.get("pgap_description", ""),
                    "tigrfam_id":          rf.get("tigrfam_id", ""),
                    "tigrfam_description": rf.get("tigrfam_description", ""),
                    "ncbifam_id":          rf.get("ncbifam_id", ""),
                    "ncbifam_description": rf.get("ncbifam_description", ""),
                    "cog_id":              rf.get("cog_id", ""),
                    "cog_description":     rf.get("cog_description", ""),
                    "pfam_id":             rf.get("pfam_id", ""),
                    "pfam_description":    rf.get("pfam_description", ""),
                    "geneprop_id":         rf.get("geneprop_id", ""),
                    "geneprop_description":rf.get("geneprop_description", ""),
                    "interpro_id":         rf.get("interpro_id", ""),
                    "interpro_description":rf.get("interpro_description", ""),
                    "kegg_id":             rf.get("kegg_id", ""),
                    "kegg_description":    rf.get("kegg_description", ""),
                    "eggnog_id":           rf.get("eggnog_id", ""),
                    "eggnog_description":  rf.get("eggnog_description", ""),
                    "uniprot_id":          rf.get("uniprot_id", ""),
                    "uniprot_description": rf.get("uniprot_description", ""),
                    # Specialized databases
                    "tcdb_id":             rf.get("tcdb_id", ""),
                    "tcdb_description":    rf.get("tcdb_description", ""),
                    "merops_id":           rf.get("merops_id", ""),
                    "merops_description":  rf.get("merops_description", ""),
                    "dbcan_id":            rf.get("dbcan_id", ""),
                    "dbcan_description":   rf.get("dbcan_description", ""),
                    # Localization
                    "signalp6_prediction": rf["signalp6_prediction"],
                    "signalp6_score":      rf["signalp6_score"],
                    "phobius_prediction":  rf["phobius_prediction"],
                    "phobius_topology":    rf["phobius_topology"],
                    "tmbed_prediction":    rf["tmbed_prediction"],
                    "psortb_prediction":   rf["psortb_prediction"],
                    "psortb_score":        rf["psortb_score"],
                    # Mechanical assessment
                    "mechanical_c1":    rf["mechanical_c1"],
                    "mechanical_c2":    rf["mechanical_c2"],
                    "mechanical_c3":    rf["mechanical_c3"],
                    "mechanical_c4":    rf["mechanical_c4"],
                    "mechanical_score": rf["mechanical_score"],
                    "mechanical_tier":  rf["mechanical_tier"],
                    "mechanical_flag":  rf["mechanical_flag"],
                    # LLM assessment
                    "neighborhood_fit":          parsed["neighborhood_fit"],
                    "neighborhood_contradiction": parsed["neighborhood_contradiction"],
                    "llm_c3":      f"{parsed['llm_c3']:.4f}" if parsed["llm_c3"] is not None else "",
                    "llm_c4":      f"{parsed['llm_c4']:.4f}" if parsed["llm_c4"] is not None else "",
                    "llm_c5":      f"{parsed['llm_c5']:.4f}" if parsed["llm_c5"] is not None else "",
                    "topology":    parsed["topology"],
                    "specialized_db_agreement": parsed["specialized_db_agreement"],
                    "operon_coherent":          parsed["operon_coherent"],
                    "formula_llm_score":         f"{formula_score:.4f}" if formula_score is not None else "",
                    "llm_assessment_score":      f"{parsed['llm_confidence']:.4f}" if parsed["llm_confidence"] is not None else "",
                    "llm_agreement_with_label":  parsed["verdict"],
                    "llm_confidence_on_label":   f"{formula_score:.4f}" if formula_score is not None else "",
                    "llm_reasoning_text":        llm_reasoning,
                })

                verdict_counts[parsed["verdict"]] = verdict_counts.get(parsed["verdict"], 0) + 1
                n_done += 1

    LOGGER.info(f"Assessed {n_done} genes → {reports_dir}/ + {summary_path}")
    LOGGER.info("  LLM verdict distribution:")
    for v, c in sorted(verdict_counts.items(), key=lambda kv: -kv[1]):
        LOGGER.info(f"    {v:15s} {c:6d} ({100.0 * c / n_done:.1f}%)")


if __name__ == "__main__":
    main()
