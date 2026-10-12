#!/usr/bin/env bash
#PBS -N diffusion_panel
#PBS -P 11004454
#PBS -q normal
#PBS -l select=1:ncpus=8:ngpus=2
#PBS -l walltime=24:00:00
#PBS -j oe
# ============================================================================
# Panel job — NSCC ASPIRE 2A. Submit from the repo root (~/DiffusionControl):
#
#   qsub -v PHASE=screen     submit_panel.sh   # stock screen + diagnosis
#   qsub -v PHASE=intervene  submit_panel.sh   # controlled arms (5 seeds)
#   qsub -v PHASE=degrade    submit_panel.sh   # MSA-depth boost (controls)
#   qsub -v PHASE=analyze    submit_panel.sh   # stats + verdicts
#   qsub -v PHASE=all        submit_panel.sh   # the four above, in order
#
# Everything is resumable: finished seeds print [skip], so a walltime kill
# or a crash means: fix the cause, resubmit the same PHASE. Nothing redone.
#
# Keep the ngpus= above in sync with reality — GPUS is auto-detected below.
# ============================================================================
set -euo pipefail

PHASE="${PHASE:-screen}"
case "$PHASE" in screen|intervene|degrade|analyze|all) ;; *)
    echo "FATAL: PHASE must be screen|intervene|degrade|analyze|all (got '$PHASE')" >&2
    exit 2 ;;
esac

SCRATCH="${SCRATCH:-/scratch/users/industry/ism/behjoo}"
ENV_PREFIX="${ENV_PREFIX:-$SCRATCH/envs/boltz_ctrl}"
REPO="${REPO:-$HOME/DiffusionControl}"

# caches on scratch (batch jobs do not read ~/.bashrc)
export PIP_CACHE_DIR="$SCRATCH/pipcache" TMPDIR="$SCRATCH/tmp" XDG_CACHE_HOME="$SCRATCH/cache"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ENV_PREFIX"
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib/python3.11/site-packages/nvidia/nvrtc/lib:${LD_LIBRARY_PATH:-}"

cd "$REPO"

# knobs — GPUS follows the actual allocation; chunking avoids OOM >= ~450 aa
GPUS="${GPUS:-$(nvidia-smi -L 2>/dev/null | grep -c '^GPU' || echo 1)}"
export GPUS
export BOLTZ_EXTRA_ARGS="${BOLTZ_EXTRA_ARGS:---max_parallel_samples 30}"
# frozen by the protocol (do NOT set): N STEPS SEEDS CTRL_SEEDS DEGRADE_DEPTH

echo "== panel job =="
echo "  job     : ${PBS_JOBID:-interactive}"
echo "  node    : $(hostname)   phase: $PHASE   gpus: $GPUS"
echo "  extra   : ${BOLTZ_EXTRA_ARGS:-<none>}"
nvidia-smi --query-gpu=index,name,memory.total --format=csv,noheader || true

echo "== freeze check =="
python3 panel/make_manifest.py --verify   # corrupt/partial deploys die here

case "$PHASE" in
    screen)    bash panel/run_panel.sh screen ;;
    intervene) bash panel/run_panel.sh intervene ;;
    degrade)   bash panel/run_panel.sh degrade ;;
    analyze)   bash panel/run_panel.sh analyze ;;
    all)       bash panel/run_panel.sh screen && bash panel/run_panel.sh intervene \
            && bash panel/run_panel.sh degrade && bash panel/run_panel.sh analyze ;;
esac

echo "== done: $PHASE =="
echo "REMINDER: scratch purges at 30 days — copy results/ and panel_msa/*.a3m"
echo "          back to the paper repo / archival storage now."
