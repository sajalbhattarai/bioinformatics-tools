#!/bin/bash
#SBATCH --job-name=c3fig_cache
#SBATCH --account=your-account
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=logs/c3-figures/cache-%j.out
#SBATCH --error=logs/c3-figures/cache-%j.err
#
# Builds the pickle caches (genes.pkl + operons.pkl) that every C3 figure script reads.
# Configure before submitting:
#   * set --account above to your SLURM allocation
#   * export these (or edit the defaults below):
#       C3_FIGDIR  path to .../scoring/c3_figures in your bioinformatics-tools checkout
#       C3_RUN     the pipeline output run (.../output/<run-id>)
# Run with:  sbatch submit_c3fig_cache.sbatch
set -euo pipefail

SCRIPTS="${C3_FIGDIR:?set C3_FIGDIR to .../scoring/c3_figures}"
RUN="${C3_RUN:?set C3_RUN to your .../output/<run-id>}"
OUT="${RUN}/scoring/c3-genes-comprehensive-stats"

mkdir -p logs/c3-figures

echo "[$(date -Is)] building cache ..."
python3 "${SCRIPTS}/c3fig_00_build_cache.py" \
    --labeling-root "${RUN}" \
    --output-dir "${OUT}"
echo "[$(date -Is)] cache build complete"
cat "${OUT}/_cache/cache-manifest.txt"
