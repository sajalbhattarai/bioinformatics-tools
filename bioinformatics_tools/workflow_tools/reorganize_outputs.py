#!/usr/bin/env python3
"""reorganize_outputs.py -- final per-organism output cleanup for a margie run.

Runs ONCE at the end of the batch (after the run-level global report), reducing
each organism folder to just the things a user cares about:

  <run>/<organism>/
    FINAL_ANNOTATION_WITH_CONFIDENCE.tsv     (promoted out of scoring/)
    FINAL_ANNOTATION_WITH_CONFIDENCE.xlsx    (colored workbook, generated here)
    <figures-dirname>/                        (was scoring/figures/, e.g. "diagrams")
    per-tool-phased-output/                   (EVERYTHING else -- all tool phase
                                               folders + the rest of scoring/)

Why it runs at the very end (not per-organism during the batch): the run-level
pangenome report (make_global_report.py / run_report_figures_global) reads EVERY
organism's scoring/scored-labeled-genes-confidence-final.tsv and
labeling/labeled-genes.tsv directly off disk. Moving those before it runs would
break it, so this step is invoked only after that global report finishes.

Idempotent and parallel: safe to re-run; organisms are processed concurrently
(the per-organism Excel generation is the slow part and is what parallelises).
This is a plain post-processing script, NOT a Snakemake rule -- moving rule
outputs would confuse Snakemake's completeness tracking on a re-run.
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

FINAL_TSV = "FINAL_ANNOTATION_WITH_CONFIDENCE.tsv"
FINAL_XLSX = "FINAL_ANNOTATION_WITH_CONFIDENCE.xlsx"
PTP = "per-tool-phased-output"
# run-level (not per-organism) folders that must never be treated as an organism
_RUN_LEVEL = {"scoring", "sqlite", "ani", "aai", "closest", "mauve",
              "original_container_outputs", "logs", "genome_pool"}


def _is_organism_dir(d: Path) -> bool:
    """A per-organism folder: has scoring/FINAL... (fresh), or FINAL... at top
    plus per-tool-phased-output/ (already reorganized)."""
    if d.name in _RUN_LEVEL:
        return False
    return ((d / "scoring" / FINAL_TSV).is_file()
            or (d / FINAL_TSV).is_file()
            or (d / PTP).is_dir())


def reorganize_one(org_dir, figures_dirname: str = "diagrams",
                   excel_script: str | None = None, python: str | None = None) -> dict:
    org = Path(org_dir)
    diag = org / figures_dirname
    ptp = org / PTP
    scoring = org / "scoring"
    result = {"organism": org.name, "actions": [], "excel": None, "error": None}
    try:
        # locate the FINAL tsv wherever it currently is (scoring/ if fresh, top if re-run)
        final_scoring = scoring / FINAL_TSV
        final_top = org / FINAL_TSV
        final_src = final_scoring if final_scoring.is_file() else (
            final_top if final_top.is_file() else None)

        # 1. colored Excel next to the FINAL tsv (skip if already there)
        if excel_script and final_src is not None:
            xlsx_dst = final_src.parent / FINAL_XLSX
            if xlsx_dst.is_file():
                result["excel"] = "exists"
            else:
                proc = subprocess.run(
                    [python or sys.executable, str(excel_script),
                     "--input", str(final_src), "--output", str(xlsx_dst)],
                    capture_output=True, text=True)
                if proc.returncode == 0:
                    result["excel"] = "generated"
                else:
                    result["excel"] = f"failed(rc={proc.returncode})"
                    result["actions"].append(
                        f"excel-gen-failed: {(proc.stderr or '').strip()[:200]}")

        # 2. scoring/figures -> <figures-dirname>/ at the organism top level
        figsrc = scoring / "figures"
        if figsrc.is_dir() and not diag.exists():
            shutil.move(str(figsrc), str(diag))
            result["actions"].append(f"figures -> {figures_dirname}/")

        # 3. promote FINAL tsv + xlsx to the top level
        for name in (FINAL_TSV, FINAL_XLSX):
            src = scoring / name
            if src.is_file() and not (org / name).exists():
                shutil.move(str(src), str(org / name))
                result["actions"].append(f"promote {name}")

        # 4. everything else -> per-tool-phased-output/
        keep = {figures_dirname, PTP, FINAL_TSV, FINAL_XLSX}
        ptp.mkdir(exist_ok=True)
        for entry in list(org.iterdir()):
            if entry.name in keep:
                continue
            dst = ptp / entry.name
            if dst.exists():
                continue  # leftover from a partial/previous run -> leave it
            shutil.move(str(entry), str(dst))
            result["actions"].append(f"-> {PTP}/{entry.name}")
    except Exception as exc:  # never raise into the pool; report per-organism
        result["error"] = f"{type(exc).__name__}: {exc}"
    return result


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-root", required=True, help="the timestamped run output dir")
    ap.add_argument("--genomes", nargs="*", default=None,
                    help="genome stems to reorganize; default = auto-detect organism dirs")
    ap.add_argument("--figures-dirname", default="diagrams",
                    help="name for the promoted scoring/figures folder (default: diagrams)")
    ap.add_argument("--excel-script", default=None,
                    help="path to make-final-excel.py; if given, a colored .xlsx is generated")
    ap.add_argument("--python", default=None, help="python for the Excel subprocess")
    ap.add_argument("--workers", type=int, default=4, help="organisms reorganized in parallel")
    args = ap.parse_args()

    run = Path(args.run_root)
    if not run.is_dir():
        print(f"[reorganize] ERROR: run root not found: {run}", file=sys.stderr)
        raise SystemExit(1)

    if args.genomes:
        org_dirs = [run / g for g in args.genomes if (run / g).is_dir()]
    else:
        org_dirs = [d for d in sorted(run.iterdir()) if d.is_dir() and _is_organism_dir(d)]

    if not org_dirs:
        print("[reorganize] no organism folders to reorganize")
        return

    print(f"[reorganize] {len(org_dirs)} organism(s); figures -> {args.figures_dirname}/; "
          f"excel={'on' if args.excel_script else 'off'}; workers={args.workers}")
    ok = 0
    with ThreadPoolExecutor(max_workers=max(1, args.workers)) as ex:
        futs = {ex.submit(reorganize_one, d, args.figures_dirname, args.excel_script, args.python): d
                for d in org_dirs}
        for fut in as_completed(futs):
            r = fut.result()
            tag = "ERROR" if r["error"] else "ok"
            print(f"  [{tag}] {r['organism']}: excel={r['excel']} "
                  f"moves={len(r['actions'])}" + (f"  ERR={r['error']}" if r["error"] else ""))
            if not r["error"]:
                ok += 1
    print(f"[reorganize] done: {ok}/{len(org_dirs)} reorganized")


if __name__ == "__main__":
    main()
