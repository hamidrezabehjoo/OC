#!/usr/bin/env bash
# One-time environment build for the panel test — NSCC ASPIRE 2A runbook.
# Repo root = this directory. Safe to re-run (every step is idempotent).
#
#   bash setup_env.sh
#
# Cluster constraints this encodes (each one cost us a failure):
#   * home is 50 GB and quota'd; scratch (100 TB) holds everything.
#     NOTE: scratch has a 30-day purge — archive results off it regularly.
#   * node driver is CUDA 12.8 -> torch MUST come from the cu128 wheel
#     index and be installed BEFORE boltz (PyPI torch is cu13 and dies
#     with "Driver too old (found version 12080)").
#   * cuequivariance needs the nvrtc lib dir on LD_LIBRARY_PATH or it
#     silently falls back to slow generic kernels.
#   * ~/.boltz checkpoints are GB-scale -> symlinked to scratch.
#   * any boltz reinstall/upgrade wipes the diffusionv2.py patch ->
#     step 7 re-applies it (this script doubles as apply_boltz_patch.sh).

set -euo pipefail

SCRATCH="${SCRATCH:-/scratch/users/industry/ism/behjoo}"
ENV_PREFIX="${ENV_PREFIX:-$SCRATCH/envs/boltz_ctrl}"
TORCH_VER="${TORCH_VER:-2.11.0}"   # cu128 wheel; node driver is CUDA 12.8
BOLTZ_VER="2.2.1"                  # pinned — the controller patch targets this

echo "== step 0: caches on scratch =="
mkdir -p "$SCRATCH"/conda/pkgs "$SCRATCH"/pipcache "$SCRATCH"/tmp "$SCRATCH"/cache "$SCRATCH"/boltz_cache
conda config --add pkgs_dirs "$SCRATCH/conda/pkgs" 2>/dev/null || true
if ! grep -q "DiffusionControl scratch caches" ~/.bashrc 2>/dev/null; then
cat >> ~/.bashrc <<EOF

# DiffusionControl scratch caches (setup_env.sh)
export PIP_CACHE_DIR="$SCRATCH/pipcache"
export TMPDIR="$SCRATCH/tmp"
export XDG_CACHE_HOME="$SCRATCH/cache"
EOF
fi
export PIP_CACHE_DIR="$SCRATCH/pipcache" TMPDIR="$SCRATCH/tmp" XDG_CACHE_HOME="$SCRATCH/cache"

echo "== step 1: conda env on scratch (never in HOME) =="
if [ -d "$ENV_PREFIX" ]; then
    echo "  exists: $ENV_PREFIX (reusing)"
else
    conda create -p "$ENV_PREFIX" python=3.11 -y
fi
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ENV_PREFIX"

echo "== step 2: torch $TORCH_VER+cu128 BEFORE boltz (non-negotiable) =="
pip install "torch==$TORCH_VER" --index-url https://download.pytorch.org/whl/cu128

echo "== step 3: boltz[cuda]==$BOLTZ_VER (torch requirement already satisfied) =="
pip install "boltz[cuda]==$BOLTZ_VER" "numpy==1.26.4" scipy

echo "== step 4: nvrtc lib on LD_LIBRARY_PATH (cue kernels) =="
NVR="$CONDA_PREFIX/lib/python3.11/site-packages/nvidia/nvrtc/lib"
export LD_LIBRARY_PATH="$NVR:${LD_LIBRARY_PATH:-}"
if ! grep -q "nvidia/nvrtc/lib" ~/.bashrc 2>/dev/null; then
    echo "export LD_LIBRARY_PATH=\"$NVR:\$LD_LIBRARY_PATH\"  # DiffusionControl cue kernels" >> ~/.bashrc
fi
echo "  (also add this export to every PBS script that runs boltz)"
python -c "import cuequivariance_ops; print('cue ops OK')"

echo "== step 5: boltz model cache on scratch =="
if [ -d "$HOME/.boltz" ] && [ ! -L "$HOME/.boltz" ]; then
    mv "$HOME/.boltz"/* "$SCRATCH/boltz_cache/" 2>/dev/null || true
    rmdir "$HOME/.boltz" 2>/dev/null || true
fi
[ -L "$HOME/.boltz" ] || ln -s "$SCRATCH/boltz_cache" "$HOME/.boltz"

echo "== step 6: torch sanity (cuda=False is normal on a login node) =="
python -c "import torch; print('cuda:', torch.cuda.is_available(), '| torch', torch.__version__, '| built cu', torch.version.cuda)"

echo "== step 7: controller patch (pinned file, md5-gated) =="
md5sum diffusionv2.py   # must print 284a8d4e58da9bee4defe972d686b22a
BOLTZ_PKG=$(python -c "import boltz, os; print(os.path.dirname(boltz.__file__))")
[ -d "$BOLTZ_PKG/model/modules" ] || { echo "FATAL: boltz package not found - env active?"; exit 1; }
[ -f "$BOLTZ_PKG/model/modules/diffusionv2.py.orig" ] || \
    cp "$BOLTZ_PKG/model/modules/diffusionv2.py" "$BOLTZ_PKG/model/modules/diffusionv2.py.orig"
[ "$(grep -c NU_ALPHA "$BOLTZ_PKG/model/modules/diffusionv2.py.orig")" = "0" ] || \
    echo "  WARNING: .orig backup is not pristine stock (was taken after an earlier patch)"
cp diffusionv2.py "$BOLTZ_PKG/model/modules/diffusionv2.py"
diff "$BOLTZ_PKG/model/modules/diffusionv2.py" diffusionv2.py && echo "PATCH APPLIED"
python -c "import torch, boltz.model.modules.diffusionv2 as d; \
print('loads from:', d.__file__); \
print('NU_ALPHA hook:', hasattr(d, '_NU_ALPHA')); \
print('SVD guard live:', torch.linalg.svd.__name__ == '_guarded_linalg_svd')"
boltz --help > /dev/null && echo "boltz CLI ok"

echo "== step 8: freeze env spec (this is what survives the 30-day purge) =="
mkdir -p "$HOME/env-specs"
pip freeze > "$HOME/env-specs/boltz_ctrl_freeze_$(date +%F).txt"
echo "  wrote $HOME/env-specs/boltz_ctrl_freeze_$(date +%F).txt"

cat <<EOF

DONE. Next:
  conda activate $ENV_PREFIX
  python3 panel/make_manifest.py --verify        # must print: manifest OK

Compute-node smoke test (never compute on the login node):
  qsub -I -l select=1:ngpus=1 -l walltime=2:00:00 -P 11004454 -q normal
  conda activate $ENV_PREFIX
  python -c "import torch; print(torch.cuda.is_available(), torch.cuda.get_device_name(0))"
  # expect: True NVIDIA A100-SXM4-40GB ... then, from the repo root:
  boltz predict yamls/adk.yaml --use_msa_server --use_potentials \
      --output_format pdb --seed 42 --out_dir results_smoke --override

Then the panel: ./panel/run_panel.sh screen
EOF
