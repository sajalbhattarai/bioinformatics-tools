#!/bin/bash
#SBATCH --job-name=c3-figures
#SBATCH --account=your-account
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:30:00
#SBATCH --output=logs/c3-figures/c3fig_all_%j.out
#SBATCH --error=logs/c3-figures/c3fig_all_%j.err
#
# Regenerate the C3 Operon Context Confidence (OCC) figure suite (cache + 7 figures + TSVs).
# Uses system python3 (numpy/pandas/scipy/matplotlib + a Times-compatible serif).
#
# Configure before submitting:
#   * set --account above to your SLURM allocation
#   * export these (or edit the defaults below):
#       C3_FIGDIR  path to .../scoring/c3_figures in your bioinformatics-tools checkout
#       C3_RUN     the pipeline output run to score (.../output/<run-id>)
# SLURM writes the logs above relative to your submit directory.

set -u

FIGDIR="${C3_FIGDIR:?set C3_FIGDIR to .../scoring/c3_figures}"
RUN="${C3_RUN:?set C3_RUN to your .../output/<run-id>}"
OUT="$RUN/scoring/c3-genes-comprehensive-stats"

mkdir -p logs/c3-figures

cd "$FIGDIR" || exit 1
echo "[c3fig] host=$(hostname) start=$(date)"
echo "[c3fig] python=$(python3 --version 2>&1)"

# 1) Rebuild the primary cache (genes / operons / adjacent pairs) from scratch.
echo "[c3fig] ===== building cache ====="
python3 c3fig_00_build_cache.py --labeling-root "$RUN" --output-dir "$OUT"
rc=$?
if [ $rc -ne 0 ]; then echo "[c3fig] CACHE BUILD FAILED (rc=$rc)"; exit $rc; fi

# 2) Build the dynamic Operon Context Confidence (OCC) reference database from the
#    cache (persisted to _cache/occ_reference.pkl). This is the production scoring
#    artifact; fig 07 reads it. To fold in a NEWLY-LABELED organism later without a
#    full rebuild, run:  python3 score_occ.py update --stats-dir "$OUT" --organisms <org>
echo "[c3fig] ===== building OCC reference database ====="
python3 score_occ.py build --stats-dir "$OUT"
rc=$?
if [ $rc -ne 0 ]; then echo "[c3fig] OCC DB BUILD FAILED (rc=$rc)"; exit $rc; fi

# 3) Run every figure script (01..07) in order. Co-occurrence / OCC-reference
#    caches are built lazily by the first figure that needs them.
fail=0
for script in $(ls c3fig_[0-9][0-9]_*.py | grep -v c3fig_00_ | sort); do
    echo "[c3fig] ===== $script ====="
    python3 "$script" --stats-dir "$OUT"
    rc=$?
    if [ $rc -ne 0 ]; then echo "[c3fig] !!! $script FAILED (rc=$rc)"; fail=$((fail+1)); fi
done

# 4) Write the per-gene OCC scores table (production scoring output).
echo "[c3fig] ===== writing per-gene OCC scores ====="
python3 score_occ.py score --stats-dir "$OUT" -o "$OUT/per_gene_occ.tsv"

npng=$(find "$OUT/figures" -name '*.png' | wc -l)
ntsv=$(find "$OUT/figures" -name '*.tsv' | wc -l)
echo "[c3fig] done=$(date) png=$npng tsv=$ntsv failed_scripts=$fail"
exit $fail
