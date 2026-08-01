#!/bin/bash
#SBATCH --job-name=margie-chat-llm
#SBATCH --account=lindems
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=%x-%j.log
#
# ---------------------------------------------------------------------------
# Long-lived inference endpoint for "Chat about the genome".
#
# Mirrors phase15 (rule run_llm) exactly: gpu partition, gpu:1, 32G, the same
# llm.sif container, and the same PYTORCH_HIP_ALLOC_CONF export -- the GPUs
# here are AMD/ROCm, which is also why this is not vLLM or TGI (neither is in
# the container, and neither is reliable on ROCm).
#
# The difference from phase15: that scores a whole genome and exits; this holds
# the model in memory and answers questions until walltime.
#
#   sbatch chat-server.sbatch
#
# On startup it writes host/port to $ADVERTISE. The API reads that file to find
# the endpoint; if the file is absent, chat reports itself offline rather than
# erroring. When this job ends the file is removed, so the UI degrades cleanly.
#
# Walltime is the real constraint: at 8h the job dies and chat goes offline
# until resubmitted. Raise --time, or add a cron/systemd-style resubmit, if it
# should stay up unattended.
# ---------------------------------------------------------------------------
set -uo pipefail

BACKEND="${MARGIE_BACKEND:-/scratch/negishi/bhattar3/margie-backend/bioinformatics-tools}"
SIF="${MARGIE_LLM_SIF:-/scratch/negishi/bhattar3/margie/sif/llm.sif}"
# The BASE model, deliberately -- not fused-model. Chat is a separate feature
# from phase15: phase15 uses the LoRA-fused model because it was fine-tuned to
# emit scores, and that is exactly the wrong bias here. This feature has to say
# "the evidence does not answer that" rather than produce a confident number,
# and refusing was never a rewarded behaviour for the fuse.
MODEL="${MARGIE_LLM_MODEL:-/scratch/negishi/bhattar3/margie/db/llm/base}"
# Must be on a shared filesystem: this job runs on a GPU node, the API runs on
# a login node. /tmp is node-local and would be invisible to the reader.
ADVERTISE="${MARGIE_CHAT_ADVERTISE:-$HOME/.local/share/bsp/chat-endpoint.json}"

echo "=== MARGIE chat server ==="
echo "node      : $(hostname -f)"
echo "job       : ${SLURM_JOB_ID:-<interactive>}"
echo "model     : $MODEL"
echo "container : $SIF"
echo "advertise : $ADVERTISE"
echo "start     : $(date)"

for p in "$SIF" "$MODEL" "$BACKEND"; do
  [ -e "$p" ] || { echo "MISSING: $p" >&2; exit 2; }
done

# ROCm allocator setting, same as phase15.
export PYTORCH_HIP_ALLOC_CONF=expandable_segments:True

# --bind the backend so the container can see chat_server.py, and the model dir
# in case it lives outside the container's default binds.
apptainer exec --nv \
    --bind "$BACKEND":"$BACKEND" \
    --bind "$(dirname "$MODEL")":"$(dirname "$MODEL")" \
    --bind "$(dirname "$ADVERTISE")":"$(dirname "$ADVERTISE")" \
    "$SIF" \
    python3 "$BACKEND/bioinformatics_tools/workflow_tools/llm/chat_server.py" \
        --model "$MODEL" \
        --advertise "$ADVERTISE"

rc=$?
echo "finish    : $(date)  rc=$rc"
# Belt and braces: chat_server removes its own advert on clean exit, but not if
# SLURM kills it at walltime. A stale advert would leave the UI pointing at a
# dead endpoint.
[ -f "$ADVERTISE" ] && python3 -c "
import json,os,sys
p=os.environ['ADVERTISE']
try:
    if json.load(open(p)).get('host','').startswith(os.uname().nodename.split('.')[0]):
        os.unlink(p); print('stale advert removed')
except Exception: pass
" 2>/dev/null
exit $rc
