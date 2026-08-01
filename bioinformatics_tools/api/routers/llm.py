"""Genome chat -- evidence-grounded Q&A over a finished run.

Separate from phase15 (workflow_tools/llm/score-genes-llm.py), which *produces*
scores as part of a run using the LoRA-fused model. This only *explains*
results that already exist, and it runs the BASE model: the fuse was trained to
emit a score, and the behaviour needed here is the opposite -- say "the evidence
does not answer that" instead of producing a confident number.

Route:  POST /v1/llm/chat  -> {"answer": str, "context_summary": {...}}
        GET  /v1/llm/status -> whether the inference endpoint is up

Grounding: the browser sends only an IDENTIFIER (job + organism + optional
gene/operon). This module rebuilds the context server-side from that job's own
FINAL table, so a tampered page cannot widen its own grounding -- constraining
the context is the entire point of the feature.

The model is reached via chat_server.py running on a GPU node, discovered
through the advert file it publishes. No advert = chat is offline, which is
reported as a clean 503 rather than an error.
"""
from __future__ import annotations

import csv
import io
import json
import logging
import os
import re
import urllib.error
import urllib.request
from pathlib import Path

from fastapi import APIRouter, Depends, HTTPException
from fastapi.responses import StreamingResponse
from pydantic import BaseModel

from bioinformatics_tools.api.auth import get_current_user
from bioinformatics_tools.api.routers.ssh import _build_connection, _resolve_job_work_dir
from bioinformatics_tools.utilities import ssh_sftp

LOGGER = logging.getLogger(__name__)

router = APIRouter(prefix="/v1/llm", tags=["llm"])

ADVERT_PATH = Path(os.getenv(
    "MARGIE_CHAT_ADVERTISE",
    os.path.expanduser("~/.local/share/bsp/chat-endpoint.json"),
))

FINAL_TSV = "FINAL_ANNOTATION_WITH_CONFIDENCE.tsv"
TIERS = ["highest", "high", "medium", "fair", "low"]

# Carried over from score-genes-llm.py's SYSTEM_PROMPT where it still applies.
# Those scoping rules encode real failure modes this pipeline already hit --
# neighbour evidence being attributed to the candidate gene, and EC granularity
# being read as conflict -- so they are repeated verbatim in intent rather than
# reinvented. What is new is the grounding contract and the citation duty.
SYSTEM_PROMPT = """\
You are explaining the results of a MARGIE(SB) genome annotation to the person \
who ran it. You are given that pipeline's own records. Explain what they mean. \
Do not re-score anything.

GROUNDING CONTRACT -- this overrides everything else:
  * Every claim must be traceable to a specific field or row in the records below.
  * You have no outside knowledge. You do not know what any gene "usually" does, \
what any organism is "known for", or what any database "normally" reports. If it \
is not in the records, you do not know it.
  * If the records do not answer the question, say what is missing and stop. Do \
not substitute plausible biology.
  * Absent is not zero, not average, and not "typical". Never infer a missing value.
  * Cite the source of each claim inline -- the column name or the tool name.

SCOPING RULES (these are real error modes, not hypotheticals):
  * The gene asked about is the CANDIDATE, and it is the one under "THE GENE \
THIS QUESTION IS ABOUT". Everything under "OTHER RECORDS" is a DIFFERENT GENE. \
Never report another gene's descriptor, scores, operon or EC as the candidate's \
-- if you mention one, name it. Answering about the wrong gene is the single \
worst failure here: it is confident, detailed, and completely wrong.
  * EC specificity differences are NOT conflicts. 2.7.1.2 vs 2.7.1.- is the same \
enzyme at different resolution. Only different reaction classes (2.7.1.2 vs \
1.1.1.27) are genuine conflicts.
  * C1 (tool coverage) and C2 (operon reliability) are arithmetic. Report them as \
given; do not re-derive them.
  * "Different function" is not "contradiction". Uninformative or hypothetical \
hits carry no signal either way -- say so rather than reading them as negative.

OPERONS: when members of an operon are supplied, you may reason about what the \
unit does -- read the members' product descriptors, order and strands together. \
That reasoning is INTERPRETATION and belongs under that heading. The members, \
their products, their coordinates and the operon id are evidence and belong \
under From the evidence. Keep the two apart: "these five members are named as \
flagellar basal-body components (evidence)" is a different kind of statement \
from "so this operon likely encodes a flagellar assembly module \
(interpretation)".

VOCABULARY:
  * Tiers, most to least confident: highest, high, medium, fair, low. \
NOT_APPLICABLE_NON_CODING marks non-coding features.
  * Two tier columns exist -- CONFIDENCE_TIER and CONFIDENCE_TIER_hybrid. State \
which one you are quoting.
  * Confidence values are decimals 0-1, never percentages.
  * NEEDS_REVIEW means the pipeline wants human attention. It does NOT mean the \
annotation is wrong.

CONVERSATION: earlier turns may be supplied so that references like "it" or \
"that gene" resolve. They are context for understanding the question ONLY. \
Never cite the conversation as evidence, and never treat something you said \
earlier as established fact -- if a claim matters, it must come from the \
records supplied now. If a follow-up refers to a gene whose record is not in \
this turn's records, say so rather than answering from memory of the last turn.

WHAT YOU ARE LOOKING AT:
  You are given (a) the COLUMN SCHEMA of this run's FINAL_ANNOTATION_WITH_
CONFIDENCE.tsv -- what each field means -- (b) genome-level aggregates, and \
(c) the FULL RECORDS of any genes relevant to the question, when the question \
names one. Study the records before answering. Aggregates describe the WHOLE \
GENOME and say nothing about any individual gene: never read a genome total as \
a statement about one gene. If a gene's record is not present, you do not know \
that gene's tier, score or review status -- say so and ask for the gene name.

ANSWER FORMAT -- every answer uses these headings, and omits any that do not apply:

  **Answer** — the direct answer, first, in one or two sentences.

  **From the evidence** — each supporting fact with the field or tool it came \
from, e.g. "CONFIDENCE_TIER = high", "NEEDS_REVIEW? = yes". Only things \
literally present in the records.

  **Interpretation** — anything you inferred, reasoned about, or drew on \
general biological knowledge to say. State plainly that it is interpretation, \
not a pipeline output. If you used no interpretation at all, omit this heading \
entirely rather than writing "none".

  **Not in the evidence** — what the records cannot answer, named specifically.

This separation is the point: the reader must be able to tell, at a glance, \
which sentences are the pipeline's findings and which are your reading of them. \
Never blend the two.

LENGTH: as long as the question genuinely needs, up to the token budget. A \
one-line question gets a few lines. Do not pad, do not restate the question, \
do not list facts that were not asked for."""


class ChatRequest(BaseModel):
    job_id: str
    organism: str                 # organism folder name inside the run
    question: str
    gene_id: str | None = None    # optional: narrows context to one gene
    max_tokens: int = 1000
    stream: bool = False
    # Recent turns, oldest first: [{"role": "you"|"margie", "text": ...}].
    # Capped server-side -- an unbounded transcript would undo the prompt
    # trimming and slow every answer down again.
    history: list[dict] = []
    # The gene the previous answer was about. Sent back by the client so a
    # follow-up ("tell me about the gene", "what does its operon say?") stays on
    # the same subject. Guessing it from the transcript is unreliable: the
    # concatenated history matches many genes and the top hit is arbitrary.
    subject_gene_id: str | None = None


def _endpoint() -> dict:
    """Current inference endpoint, or 503 if chat is not running."""
    if not ADVERT_PATH.is_file():
        raise HTTPException(
            status_code=503,
            detail="Genome chat is offline — no inference server is running. "
                   "Start it with: sbatch chat-server.sbatch",
        )
    try:
        adv = json.loads(ADVERT_PATH.read_text())
        return {"url": f"http://{adv['host']}:{adv['port']}", "model": adv.get("model", "")}
    except Exception as exc:
        raise HTTPException(status_code=503, detail=f"Chat endpoint unreadable: {exc}")


def _read_final(job_id: str, organism: str, current_user: dict) -> list[dict]:
    """Read this organism's FINAL table off the user's cluster."""
    conn = _build_connection(current_user)
    work_dir = _resolve_job_work_dir(job_id, current_user, conn)
    # Post-reorganize the FINAL table sits at the organism top level; mid-run it
    # is still under scoring/. Try both rather than assuming a layout.
    for rel in (f"{organism}/{FINAL_TSV}", f"{organism}/scoring/{FINAL_TSV}"):
        try:
            raw = b"".join(ssh_sftp.stream_remote_file(f"{work_dir}/{rel}", connection=conn))
            return list(csv.DictReader(io.StringIO(raw.decode("utf-8", "replace")),
                                       delimiter="\t"))
        except FileNotFoundError:
            continue
        except Exception as exc:
            raise HTTPException(status_code=500, detail=f"Could not read {rel}: {exc}")
    raise HTTPException(status_code=404,
                        detail=f"No {FINAL_TSV} for organism '{organism}' in this job.")


def _col(row: dict, name: str) -> str:
    """FINAL tables carry prefixed headers ('Column-AS: CONFIDENCE_TIER')."""
    if name in row:
        return (row[name] or "").strip()
    for k, v in row.items():
        if k and k.split(":")[-1].strip() == name:
            return (v or "").strip()
    return ""


# Fields worth naming explicitly. The model cannot reason about a table whose
# columns it has never been told the meaning of -- without this it treats every
# header as opaque and falls back on genome totals.
_SCHEMA_HELP = {
    "gene_id": "gene identifier, <accession>_<start><strand><len>",
    "organism_name": "genome this gene belongs to",
    "best_consensus_product_descriptor": "the product name the pipeline settled on",
    "best_consensus_product_descriptor_source": "which database that name came from",
    "BEST_PRODUCT_DESCRIPTOR(copied_here_for_convenience)": "same descriptor, duplicated for convenience",
    "RAST_feature_id": "RAST feature identifier",
    "FEATURE_TYPE": "CDS, RNA, etc.",
    "EC_EVIDENCE_STATUS": "state of EC-number evidence for this gene",
    "RAST_start": "start coordinate on its replicon",
    "RAST_end": "end coordinate",
    "RAST_strand": "+ or -",
    "CONFIDENCE_TIER": "highest | high | medium | fair | low | NOT_APPLICABLE_NON_CODING",
    "CONFIDENCE_TIER_hybrid": "same scale, operon-context-adjusted variant",
    "NEEDS_REVIEW?": "yes/no — pipeline wants a human to look",
    "IS_IN_OPERON?": "yes/no",
    "UniOP_OPERON_id": "operon identifier, or a NOT_* sentinel",
    "C1_score_tool_coverage": "fraction of tools returning an informative hit",
    "C2_score_operon_probability": "geometric mean of operon-member probabilities",
    "C3_score_operon_context": "operon/neighbourhood coherence",
    "C4_score_EC_conflict": "EC-number conflict component",
    "confidence_score": "final combined confidence, 0-1",
}


def _schema_block(rows: list[dict]) -> str:
    """Tell the model what this run's FINAL table actually contains."""
    if not rows:
        return ""
    names = []
    for k in rows[0].keys():
        if not k:
            continue
        bare = k.split(":")[-1].strip()
        names.append(f"  {bare}" + (f" — {_SCHEMA_HELP[bare]}" if bare in _SCHEMA_HELP else ""))
    return ("FINAL_ANNOTATION_WITH_CONFIDENCE.tsv — one row per gene, "
            f"{len(rows)} rows, {len(names)} columns:\n" + "\n".join(names))


# Words that must never drive retrieval. Two groups:
#   * ordinary question words
#   * SCHEMA vocabulary -- "operon" appears in every operonic gene's
#     UniOP_OPERON_id, "descriptor"/"tier"/"confidence" name columns. Matching on
#     these selects an essentially arbitrary gene, which is exactly how "what
#     does operon say about this gene?" landed on an unrelated formate
#     dehydrogenase instead of the gene under discussion.
_STOPWORDS = {
    "what", "which", "does", "have", "this", "that", "gene", "genes", "genome",
    "about", "there", "any", "can", "you", "tell", "give", "more", "into",
    "details", "detail", "look", "the", "and", "for", "with", "from", "its",
    "it", "is", "are", "has", "how", "many", "much", "please", "called",
    "then", "okay", "also", "same", "see", "mean", "say", "says", "know",
    # schema / domain vocabulary
    "operon", "operons", "operonic", "descriptor", "descriptors", "tier",
    "tiers", "confidence", "score", "scores", "review", "reviewed", "flag",
    "flagged", "product", "products", "coding", "non", "protein", "proteins",
    "annotation", "annotated", "evidence", "database", "databases", "tool",
    "tools", "hit", "hits", "value", "values", "column", "columns", "final",
}
# A token this short can substring-match almost anything -- a typo like "abot"
# matched an unrelated descriptor and made the code believe it had resolved the
# gene, which suppressed the conversation fallback entirely. Retrieval only
# trusts tokens long enough to be distinctive.
_MIN_TOKEN = 5


def _search_genes(rows: list[dict], question: str, limit: int = 6) -> list[dict]:
    """Rows whose label/id matches meaningful words in the question.

    Without this the model is handed genome totals and nothing else, and when
    asked about a named gene it invents a per-gene answer from an aggregate --
    the exact failure this is here to remove.
    """
    q = question.lower()
    words = {w for w in re.findall(r"[A-Za-z0-9_.\-]+", q)
             if len(w) >= _MIN_TOKEN and w not in _STOPWORDS}

    # Any handle should reach the same record, so coordinates count as a handle
    # too: "the gene at 1,234,500" or "genes between 10000 and 20000". Commas
    # are stripped first -- people write positions with thousands separators.
    coords = [int(n) for n in re.findall(r"\d{3,}", q.replace(",", ""))]

    scored = []
    for r in rows:
        # Every textual handle the row carries: id, label, operon, feature id,
        # EC numbers, locus tag. Matching on one of these must land the same
        # row as matching on any other.
        hay = " ".join(filter(None, (
            _col(r, "gene_id"),
            # The product name really is best_consensus_product_descriptor --
            # there is no canonical_label column in this table. Column-AW is
            # the same value copied for convenience; both are matched so either
            # spelling of the question lands the row.
            _col(r, "best_consensus_product_descriptor"),
            _col(r, "BEST_PRODUCT_DESCRIPTOR(copied_here_for_convenience)"),
            _col(r, "UniOP_OPERON_id"), _col(r, "RAST_feature_id"),
            _col(r, "EC_EVIDENCE_STATUS"),
        ))).lower()
        # Weight by token length: "glycosyl_hydrolase_malt_phosph" is far
        # stronger evidence than a 5-letter fragment.
        hits = sum(len(w) for w in words if w in hay) if hay else 0

        # Positional match: a coordinate falling inside the gene is a strong
        # signal, so it outranks a loose word match.
        if coords:
            try:
                s, e = int(_col(r, "RAST_start")), int(_col(r, "RAST_end"))
                lo, hi = min(s, e), max(s, e)
                if any(lo <= c <= hi for c in coords):
                    hits += 40   # a coordinate inside the gene is decisive
            except (ValueError, TypeError):
                pass

        if hits:
            scored.append((hits, r))
    scored.sort(key=lambda x: -x[0])
    return [r for _, r in scored[:limit]]


# Context-only records exist so the model can tell "that is a different gene".
# They do not need 49 columns each to do that job -- rendering them in full was
# 62k of a 65k-character prompt, i.e. the entire latency problem.
_BRIEF_FIELDS = ("gene_id", "best_consensus_product_descriptor", "CONFIDENCE_TIER",
                 "NEEDS_REVIEW?", "UniOP_OPERON_id", "RAST_start", "RAST_end")


def _render_brief(r: dict) -> str:
    bits = []
    for f in _BRIEF_FIELDS:
        v = _col(r, f)
        if v:
            bits.append(f"{f}={v}")
    return "    " + "; ".join(bits)


def _render_gene(r: dict) -> str:
    out = []
    for k, v in r.items():
        if not k:
            continue
        v = (v or "").strip()
        if not v:
            continue
        bare = k.split(":")[-1].strip()
        if bare == "RAST_na_sequence":
            v = f"<{len(v)} nt, omitted>"
        elif len(v) > 600:
            # Audit trails and reasoning strings run to thousands of characters.
            # Keep the head -- that is where the substance is -- and say it was cut
            # so the model does not treat the tail as absent.
            v = v[:600] + f" …[truncated, {len(v)} chars total]"
        out.append(f"    {bare}: {v}")
    return "\n".join(out)


def _genome_context(rows: list[dict], organism: str) -> tuple[str, dict]:
    """Whole-genome scope: aggregates only. Per-gene records cannot fit."""
    tally = {t: 0 for t in TIERS}
    noncoding = flagged = operonic = 0
    for r in rows:
        t = _col(r, "CONFIDENCE_TIER")
        if t in tally:
            tally[t] += 1
        elif t.startswith("NOT_"):
            noncoding += 1
        if _col(r, "NEEDS_REVIEW?").lower() == "yes":
            flagged += 1
        if _col(r, "IS_IN_OPERON?").lower() == "yes":
            operonic += 1
    lines = [f"GENOME: {organism}", f"total genes: {len(rows)}",
             "confidence tier counts (CONFIDENCE_TIER):"]
    lines += [f"  {t}: {tally[t]}" for t in TIERS]
    lines += [f"  NOT_APPLICABLE_NON_CODING: {noncoding}",
              f"genes flagged NEEDS_REVIEW: {flagged}",
              f"genes in an operon (IS_IN_OPERON?): {operonic}"]
    summary = {"scope": "genome", "genes": len(rows), "flagged": flagged,
               "tiers": tally, "non_coding": noncoding}
    return "\n".join(lines), summary


def _gene_context(rows: list[dict], gene_id: str, organism: str) -> tuple[str, dict]:
    """Single-gene scope: that gene's full FINAL record, verbatim."""
    match = next((r for r in rows if _col(r, "gene_id") == gene_id), None)
    if match is None:
        raise HTTPException(status_code=404, detail=f"gene_id '{gene_id}' not in this organism.")
    lines = [f"GENOME: {organism}", f"GENE: {gene_id}", "",
             "This gene's complete record from FINAL_ANNOTATION_WITH_CONFIDENCE.tsv.",
             "Field names are given exactly as the pipeline wrote them:", ""]
    for k, v in match.items():
        if k is None:
            continue
        v = (v or "").strip()
        if not v:
            continue
        # The nucleotide sequence is long and never the subject of a question.
        if k.split(":")[-1].strip() == "RAST_na_sequence":
            v = f"<{len(v)} nt, omitted>"
        lines.append(f"  {k.split(':')[-1].strip()}: {v}")
    return "\n".join(lines), {"scope": "gene", "gene_id": gene_id,
                              "fields": sum(1 for v in match.values() if v)}


def _resolve_model_path(current_user: dict) -> str:
    """<db_root>/llm/base, where db_root is the database folder the user set in
    the GUI. Never hardcoded: the whole point of the Profile setting is that the
    database lives wherever the user put it.

    The BASE model deliberately, not fused-model — see this module's docstring.
    """
    conn = _build_connection(current_user)
    home = current_user["home_dir"]
    cfg = {}
    try:
        cfg = ssh_sftp.read_remote_yaml(
            f"{home}/.config/bioinformatics-tools/config.yaml", connection=conn) or {}
    except Exception as exc:
        LOGGER.warning("Could not read config for db_root: %s", exc)

    # Same precedence the workflow uses: an explicit db.llm wins, else
    # <workflow>.db_root/llm. db_path() itself is Snakemake-side, so the two
    # keys are read directly here rather than importing the workflow config.
    explicit = (cfg.get("db") or {}).get("llm")
    root = None
    if explicit:
        root = str(explicit)
    else:
        for wf in ("margie_sb", "margie"):
            r = (cfg.get(wf) or {}).get("db_root")
            if r:
                root = f"{str(r).rstrip('/')}/llm"
                break
    if not root:
        root = (current_user.get("db_root") or "").rstrip("/")
        root = f"{root}/llm" if root else ""
    if not root:
        raise HTTPException(
            status_code=400,
            detail="No database folder configured. Set it in Profile → database "
                   "path (the folder containing llm/base) before starting chat.",
        )
    return f"{root.rstrip('/')}/base"


@router.post("/start")
def start_chat(current_user: dict = Depends(get_current_user)):
    """Submit the GPU job that hosts the chat model — the 'Start chat' button.

    Chat needs a resident model, and a SLURM allocation is the only way to hold
    one on this cluster. Returns immediately with the job id; the model takes a
    minute or two to load, during which /status reports offline.
    """
    if ADVERT_PATH.is_file():
        st = chat_status(current_user)
        if st.get("online"):
            return {"started": False, "already_running": True, **st}

    model = _resolve_model_path(current_user)
    conn = _build_connection(current_user)
    ssh = conn.connect()
    try:
        _, out, err = ssh.exec_command(
            f"test -d {model!r} && echo OK || echo MISSING")
        if out.read().decode().strip() != "OK":
            raise HTTPException(
                status_code=400,
                detail=f"Base model not found at {model}. Check the database "
                       "folder set in Profile.",
            )
        # Backend location is wherever this package lives on the cluster.
        backend = str(Path(__file__).resolve().parents[3])
        sbatch = f"{backend}/bioinformatics_tools/workflow_tools/llm/chat-server.sbatch"
        cmd = (f"MARGIE_LLM_MODEL={model!r} MARGIE_BACKEND={backend!r} "
               f"sbatch --parsable {sbatch!r}")
        _, out, err = ssh.exec_command(cmd)
        job = out.read().decode().strip()
        stderr = err.read().decode().strip()
        if not job.isdigit():
            raise HTTPException(status_code=500,
                                detail=f"sbatch failed: {stderr or job or 'no job id'}")
    finally:
        pass  # pooled client: closing it would break concurrent requests (see SSHConnection pool)

    LOGGER.info("chat server submitted as job %s (model %s)", job, model)
    return {"started": True, "job_id": job, "model": model,
            "detail": "Chat server queued. It reports online once the model has loaded "
                      "(usually 1–3 minutes, longer if the GPU queue is busy)."}


@router.post("/stop")
def stop_chat(current_user: dict = Depends(get_current_user)):
    """End interactive mode — called when the map page is closed.

    Asks the server to exit (which frees the GPU and removes its advert). This
    is the fast path only: it never arrives if the browser crashed, the laptop
    slept, or the network dropped, so chat_server also exits on its own idle
    timeout. Never rely on this alone to release a GPU.
    """
    if not ADVERT_PATH.is_file():
        return {"stopped": False, "detail": "chat was not running"}
    try:
        adv = json.loads(ADVERT_PATH.read_text())
        req = urllib.request.Request(
            f"http://{adv['host']}:{adv['port']}/shutdown", data=b"{}",
            headers={"Content-Type": "application/json"})
        with urllib.request.urlopen(req, timeout=10):
            pass
        LOGGER.info("chat server on %s asked to stop", adv.get("host"))
        return {"stopped": True}
    except Exception as exc:
        # Unreachable means it is already gone; clear the stale advert so the
        # UI stops pointing at a dead endpoint.
        try:
            ADVERT_PATH.unlink()
        except Exception:
            pass
        return {"stopped": True, "detail": f"endpoint already gone ({exc}); advert cleared"}


@router.get("/status")
def chat_status(current_user: dict = Depends(get_current_user)):
    """Whether the chat backend is reachable — lets the UI show an honest state."""
    if not ADVERT_PATH.is_file():
        return {"online": False, "detail": "no inference server running"}
    try:
        adv = json.loads(ADVERT_PATH.read_text())
        url = f"http://{adv['host']}:{adv['port']}/health"
        with urllib.request.urlopen(url, timeout=5) as r:
            return {"online": True, "model": json.loads(r.read()).get("model", ""),
                    "host": adv["host"]}
    except Exception as exc:
        # Advert present but unreachable => the job died without cleaning up.
        return {"online": False, "detail": f"endpoint not responding: {exc}"}


@router.post("/chat")
def chat(body: ChatRequest, current_user: dict = Depends(get_current_user)):
    question = (body.question or "").strip()
    if not question:
        raise HTTPException(status_code=400, detail="question is required")

    ep = _endpoint()
    rows = _read_final(body.job_id, body.organism, current_user)

    # Keep the last few turns, truncated. Enough for "it"/"that gene" to resolve,
    # not enough to reintroduce the prompt bloat just removed.
    hist = [h for h in (body.history or [])
            if isinstance(h, dict) and (h.get("text") or "").strip()][-6:]
    hist_text = "\n".join(
        f"{'USER' if h.get('role') == 'you' else 'YOU'}: {(h.get('text') or '')[:800]}"
        for h in hist)
    # What the user has said, for resolving references in the current question.
    prior_user = " ".join((h.get("text") or "") for h in hist
                          if h.get("role") == "you")

    if body.gene_id:
        context, summary = _gene_context(rows, body.gene_id, body.organism)
        matched = []
    else:
        context, summary = _genome_context(rows, body.organism)
        # Pull in the full record of any gene the question names. Aggregates
        # alone caused the model to answer per-gene questions from genome
        # totals ("the gene is not flagged" from a count of 1212 flagged).
        matched = _search_genes(rows, question)
        # A follow-up ("look more into it", "what about its neighbours") names
        # no gene, so searching the question alone finds nothing and the answer
        # silently falls back to genome aggregates -- which is how the earlier
        # wrong-gene answers happened. Retry against what the user asked before.
        if not matched and body.subject_gene_id:
            # Most reliable: the client tells us what the last answer was about.
            keep = [r for r in rows if _col(r, "gene_id") == body.subject_gene_id]
            if keep:
                matched = keep
                LOGGER.info("follow-up stayed on subject gene %s", body.subject_gene_id)
        if not matched and prior_user:
            # Fallback for a client that does not track the subject, or a first
            # follow-up after a page reload.
            matched = _search_genes(rows, prior_user)
            if matched:
                LOGGER.info("resolved follow-up via prior turns -> %s",
                            _col(matched[0], "gene_id"))
        # An operon can only be interpreted from its MEMBERS -- their products,
        # order and strands are the whole signal. If a matched gene sits in an
        # operon, pull in its siblings so the model reasons about the unit
        # rather than one gene in isolation.
        if matched:
            member_of = {_col(r, "UniOP_OPERON_id") for r in matched}
            member_of = {o for o in member_of if o.startswith("operon_")}
            if member_of:
                have = {id(r) for r in matched}
                siblings = [r for r in rows
                            if _col(r, "UniOP_OPERON_id") in member_of
                            and id(r) not in have]
                # Cap: a large operon would otherwise crowd out the question.
                matched = matched + siblings[:24]
        if matched:
            # The first match is the best-scoring one. Label it as THE subject
            # and everything else as background, or the model blends several
            # genes into a single answer -- observed in use: asked about DNA
            # helicase RecQ, it described a pilus assembly protein from another
            # operon, complete with that gene's scores.
            primary, others = matched[0], matched[1:]
            context += (
                "\n\nTHE GENE THIS QUESTION IS ABOUT — answer about THIS record "
                "and no other:\n"
                f"\n  [{_col(primary, 'gene_id') or '?'}] "
                f"{_col(primary, 'best_consensus_product_descriptor')}\n"
                f"{_render_gene(primary)}\n")
            if others:
                context += (
                    "\n\nOTHER RECORDS, for context only. These are DIFFERENT "
                    "GENES. Never report their values as the subject gene's, and "
                    "name the gene explicitly whenever you mention one:\n")
                for r in others:
                    context += f"\n{_render_brief(r)}\n"
            summary["matched_genes"] = [_col(r, "gene_id") for r in matched]
            # The client echoes this back as subject_gene_id on the next turn.
            summary["subject_gene_id"] = _col(primary, "gene_id")
        else:
            context += ("\n\nNo gene in this genome matched the wording of the "
                        "question, so no per-gene record is included. Only the "
                        "genome-level aggregates above are available.")

    context = _schema_block(rows) + "\n\n" + context

    convo = (f"\n\n----\nCONVERSATION SO FAR (context for resolving references "
             f"like \"it\" or \"that gene\"; NOT evidence -- never cite it as a "
             f"pipeline finding):\n{hist_text}\n" if hist_text else "")

    prompt = (f"{context}{convo}\n\n"
              f"----\nQUESTION: {question}\n\n"
              "Study the records above, then answer using the required headings. "
              "Cite the field or tool behind every fact, and keep anything you "
              "inferred under Interpretation so it cannot be mistaken for a "
              "pipeline output.")

    payload = json.dumps({"system": SYSTEM_PROMPT, "prompt": prompt,
                          "max_tokens": body.max_tokens,
                          "stream": bool(body.stream)}).encode()
    req = urllib.request.Request(f"{ep['url']}/chat", data=payload,
                                 headers={"Content-Type": "application/json"})

    if not body.stream:
        try:
            with urllib.request.urlopen(req, timeout=300) as r:
                answer = json.loads(r.read()).get("text", "")
        except urllib.error.URLError as exc:
            raise HTTPException(
                status_code=503,
                detail=f"Chat backend unreachable: {exc}. It may have hit walltime.")
        except Exception as exc:
            raise HTTPException(status_code=500, detail=f"Chat failed: {exc}")
        return {"answer": answer, "context_summary": summary, "model": ep["model"]}

    # Streaming. The connection is opened here rather than inside the generator
    # so an unreachable backend still becomes a clean 503 -- once StreamingResponse
    # has begun, the status is already sent and errors can only go in the body.
    try:
        upstream = urllib.request.urlopen(req, timeout=300)
    except urllib.error.URLError as exc:
        raise HTTPException(
            status_code=503,
            detail=f"Chat backend unreachable: {exc}. It may have hit walltime.")

    def relay():
        try:
            while True:
                chunk = upstream.read(64)
                if not chunk:
                    break
                yield chunk
        except Exception as exc:
            LOGGER.warning("stream relay ended: %s", exc)
            yield f"\n[stream interrupted: {exc}]".encode()
        finally:
            upstream.close()

    return StreamingResponse(
        relay(), media_type="text/plain; charset=utf-8",
        headers={
            # The summary cannot go in the body -- that is the answer text --
            # so it rides along as a header the client reads before streaming.
            "X-Context-Summary": json.dumps(summary),
            "X-Model": ep["model"],
            "Cache-Control": "no-cache",
            "X-Accel-Buffering": "no",
        },
    )
