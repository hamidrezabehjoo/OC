#!/usr/bin/env bash
# Panel driver for the prospective screen (OC23 + MS15 [+ optional TP16]).
# Protocol: panel/prospective_test_panel_protocol.md (v2.5).
# Operator's manual: README.md at the repo root.
#
#   ./panel/run_panel.sh fetch                        # sequences, refs, yamls
#   ./panel/run_panel.sh screen                       # stock screen: 100 x 3 seeds per target
#   ./panel/run_panel.sh intervene                    # frozen arms on diagnosed targets +
#                                                     # contra-diagnosed arms + healthy controls
#   ./panel/run_panel.sh degrade                      # MSA-depth-engineered boost (all healthy controls)
#   ./panel/run_panel.sh analyze                      # panel-level summary
#
# Filters (any phase):   --only name1,name2   --set oc23|ms15|tp16
# Execution knobs (env): GPUS=1 (screen: split targets across cards)
#                        BOLTZ_EXTRA_ARGS="" (e.g. "--max_parallel_samples 30")
# Frozen (protocol; do NOT change): N=100 STEPS=200 SEEDS="42 1 2"
#                        CTRL_SEEDS="42 1 2 7 13" N_HEALTHY_CONTROLS=5
#
# Per target <name>:
#   yamls/<name>.yaml                        fetch (construct-level input)
#   panel_data/<name>_ref1.pdb, _ref2.pdb    fetch (renumbered, chain A)
#   panel_msa/<name>.a3m                     screen (MSA frozen for all arms)
#   results_<name>_{stock,ctrl,contra}_s<seed>/   boltz output
#   eval_<name>_{stock,ctrl,contra}_s<seed>.csv   eval_panel output
#   diagnosis_<name>.txt                     frozen decision table
#   panel/panel_diagnosis.csv, panel/panel_summary.csv
#
# Prerequisites: patched Boltz-2 env (diffusionv2.py must contain NU_ALPHA;
# checked before any controlled run). Frozen arms: damp NU_ALPHA=0.5,
# boost NU_ALPHA=2.5, both NU_TW=0.5. Do NOT edit arms per target.

set -uo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

PHASE="${1:-}"
ONLY=""; SET=""
shift || true
while [ $# -gt 0 ]; do
    case "$1" in
        --only) ONLY="$2"; shift 2 ;;
        --set)  SET="$2";  shift 2 ;;
        *) echo "unknown option $1" >&2; exit 1 ;;
    esac
done

N=${N:-100}; STEPS=${STEPS:-200}
SEEDS=${SEEDS:-"42 1 2"}                    # stock screen (diagnosis)
# Amendment 2026-10-05(e): controlled/contra/degrade arms run 5 seeds.
# At 3 seeds the seed-level bootstrap is degenerate (2.5th percentile =
# smallest seed shift), so "CI excludes 0" == "3/3 positive", null rate
# 1/2^3 = 12.5%. At 5 seeds the criterion is 5/5 direction + bootstrap
# CI, null 1/2^5 = 3.1%.
CTRL_SEEDS=${CTRL_SEEDS:-"42 1 2 7 13"}
N_HEALTHY_CONTROLS=${N_HEALTHY_CONTROLS:-5}   # amendment 2026-10-05b: cap 5
DEGRADE_DEPTH=32   # frozen 2026-10-05c; do not tune per target
# Execution-only knob (2026-10-06): extra flags passed verbatim to every
# boltz predict call, e.g. BOLTZ_EXTRA_ARGS="--max_parallel_samples 30"
# for the large targets that OOM a 40GB card at 100 parallel samples.
# Verified in boltz 2.x source (diffusionv2.py): init noise, per-step eps
# and random augmentation are all drawn for the FULL multiplicity outside
# the sample-chunk loop; chunking only splits the deterministic denoiser
# forward, so per-sample noise streams are unchanged. No protocol impact.
# WARNING: with 100 samples, m in {5,10,20,25,50} is a silent no-op
# (chunk count = 100 % m + 1); use 30 (->10/chunk) or 34 (->4/chunk).
BOLTZ_EXTRA_ARGS=${BOLTZ_EXTRA_ARGS:-""}

BIGCACHE="${BIGCACHE:-$HOME}"
export BOLTZ_CACHE="${BOLTZ_CACHE:-$BIGCACHE/boltz_cache}"
export TRITON_CACHE_DIR="${TRITON_CACHE_DIR:-$BIGCACHE/triton_cache}"
export TORCHINDUCTOR_CACHE_DIR="${TORCHINDUCTOR_CACHE_DIR:-$BIGCACHE/torchinductor_cache}"
export CUDA_CACHE_PATH="${CUDA_CACHE_PATH:-$BIGCACHE/cuda_cache}"
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-$BIGCACHE/xdg_cache}"
mkdir -p "$BOLTZ_CACHE" panel_msa

TARGETS=$(python3 - "$ONLY" "$SET" <<'EOF'
import csv, sys
only = set(sys.argv[1].split(",")) if sys.argv[1] else None
wantset = sys.argv[2] or None
for r in csv.DictReader(open("panel/targets.csv")):
    if r["status"] != "ready":
        continue
    if r["optional"] == "1" and wantset != "tp16":
        continue
    if only and r["name"] not in only:
        continue
    if wantset and r["set"] != wantset:
        continue
    print(r["name"], r["foldswitch"], r["expected_state"])
EOF
)
[ -n "$TARGETS" ] || { echo "no targets selected" >&2; exit 1; }

yaml_of()    { echo "yamls/$1.yaml"; }
pred_dir()   { echo "$1/boltz_results_$2/predictions/$2"; }
# boltz names output dirs after the YAML filename STEM (boltz_results_<stem>/
# predictions/<stem>): yamls/adk.msa.yaml -> boltz_results_adk.msa/... .
# Eval must follow the stem of the yaml actually fed to boltz, not the bare
# target name (bug 2026-10-08: intervene/degrade evals all failed on this).
yaml_stem() { local b; b=$(basename "$1"); echo "${b%.yaml}"; }
basin_rule() { [ "$2" = "1" ] && echo tm || echo rmsd; }

run_eval() { # $1 results dir, $2 name, $3 foldswitch, $4 tag, $5 yaml stem
    python3 panel/eval_panel.py --pred_dir "$(pred_dir "$1" "$5")" \
        --ref1 "panel_data/$2_ref1.pdb" --ref2 "panel_data/$2_ref2.pdb" \
        --basin_rule "$(basin_rule "$2" "$3")" --escape_rmsd auto \
        --tag "$4" --out "eval_$2_$4.csv"
}

patch_gate() {
    local f n
    f=$(python3 -c "import boltz, os; print(os.path.join(os.path.dirname(boltz.__file__),'model','modules','diffusionv2.py'))")
    n=$(grep -c NU_ALPHA "$f" || true)
    if [ "$n" -eq 0 ]; then
        echo "FATAL: $f has no NU_ALPHA hook - controlled arm would silently be stock." >&2
        echo "Apply the patch first:  cp diffusionv2.py $f" >&2
        exit 1
    fi
}

# freeze the MSA from the first stock run so every later arm of this target
# uses byte-identical evolutionary input (amendment 2026-10-05)
freeze_msa() { # $1 name, $2 first results dir
    [ -s "panel_msa/$1.a3m" ] && return 0
    local a3m
    a3m=$(find "$2" -name "*.a3m" | head -1 || true)
    if [ -n "$a3m" ]; then
        cp "$a3m" "panel_msa/$1.a3m"
        # boltz's a3m writer appends one trailing NUL byte (2026-10-08:
        # parse_a3m crashes on it for MSAs smaller than max_seqs). Strip
        # exactly that artifact, then verify the copy byte-for-byte.
        python3 - "$1" <<'EOF'
import sys
p = f"panel_msa/{sys.argv[1]}.a3m"
d = open(p, "rb").read()
if d.endswith(b"\x00") and d.count(b"\x00") == 1:
    open(p, "wb").write(d[:-1])
    print(f"  freeze_msa: stripped trailing NUL (boltz writer artifact) from {p}")
elif b"\x00" in d:
    print(f"  WARNING: {p} contains NUL bytes beyond the trailing artifact -", file=sys.stderr)
    print(f"  inspect before running arms of this target", file=sys.stderr)
EOF
        cmp -s "$a3m" "panel_msa/$1.a3m" || \
            echo "  NOTE: panel_msa/$1.a3m differs from source by the stripped NUL" 
        echo "  MSA frozen: $(md5sum "panel_msa/$1.a3m" | cut -c1-12) panel_msa/$1.a3m"
    else
        echo "  WARNING: no .a3m found under $2 - ctrl arms will call the MSA server again" >&2
    fi
}

# yaml with the frozen MSA wired in (falls back to plain yaml + msa server).
# Sets MSA_YAML="" when the frozen a3m is corrupt: caller must skip the
# target (never silently fall back to the server - that would break the
# no-cross-arm-MSA-drift freeze).
msa_yaml() { # $1 name -> sets globals MSA_YAML and MSA_FLAG. Must NOT be
             # called via $(...): command substitution is a subshell and
             # the flag assignment would be lost (set -u crash, 2026-10-06)
    if [ -s "panel_msa/$1.a3m" ]; then
        if python3 -c "import sys; sys.exit(1 if b'\x00' in open('panel_msa/$1.a3m','rb').read() else 0)"; then
            sed "/      sequence:/a\\      msa: panel_msa/$1.a3m" "$(yaml_of "$1")" > "yamls/$1.msa.yaml"
            MSA_YAML="yamls/$1.msa.yaml"; MSA_FLAG=""
        else
            echo "  CORRUPT: panel_msa/$1.a3m has NUL bytes - repair or re-fetch, then rerun" >&2
            MSA_YAML=""; MSA_FLAG=""
        fi
    else
        MSA_YAML=$(yaml_of "$1"); MSA_FLAG="--use_msa_server"
    fi
}

run_arm() { # $1 name, $2 foldswitch, $3 tag(ctrl/contra), $4 alpha
    local name=$1 fs=$2 tag=$3 alpha=$4
    local y mf stem
    msa_yaml "$name"; y=$MSA_YAML; mf=$MSA_FLAG; stem=$(yaml_stem "$y")
    [ -n "$y" ] || { echo "[skip] $name $tag: frozen MSA corrupt (see above)"; return 1; }
    for s in $CTRL_SEEDS; do
        if [ -s "eval_${name}_${tag}_s${s}.csv" ]; then
            echo "[skip] $name $tag seed $s: exists"; continue
        fi
        local res="results_${name}_${tag}_s${s}"
        if [ -d "$(pred_dir "$res" "$stem")" ]; then
            # predictions already on disk (e.g. eval-path bug 2026-10-08):
            # eval is deterministic given structures - skip the GPU work
            echo "[reuse] $name $tag seed $s: predictions exist, eval only"
        else
            echo "=== $name $tag alpha=$alpha tw=0.5 seed $s ==="
            NU_ALPHA="$alpha" NU_TW=0.5 \
            boltz predict "$y" $mf --use_potentials < /dev/null \
                --diffusion_samples "$N" --sampling_steps "$STEPS" $BOLTZ_EXTRA_ARGS \
                --output_format pdb --seed "$s" --cache "$BOLTZ_CACHE" \
                --out_dir "$res" --override \
                || { echo "[FAIL boltz] $name $tag seed $s"; continue; }
        fi
        run_eval "$res" "$name" "$fs" "${tag}_s${s}" "$stem" \
            || echo "[FAIL eval] $name $tag seed $s"
    done
}

case "$PHASE" in

fetch)
    python3 panel/fetch_panel.py --targets panel/targets.csv --root . ${ONLY:+--only $ONLY}
    echo "(for TP16 rows run directly: python3 panel/fetch_panel.py --targets panel/targets.csv --root . --include-optional)"
    ;;

screen)
    # Multi-GPU execution (2026-10-06, execution-only): GPUS=N splits the
    # selected targets across N cards, one boltz process per card
    # (CUDA_VISIBLE_DEVICES pinned per worker). boltz never pools memory
    # across GPUs, so large targets still need BOLTZ_EXTRA_ARGS; this only
    # buys wall-clock. Targets are dealt greedily by length (heaviest to
    # the lightest worker) so the 700-850 aa proteins don't share a card.
    N_GPUS=${GPUS:-1}
    TGT="$TARGETS" python3 - "$N_GPUS" <<'EOF'
import csv, os, sys
n = int(sys.argv[1])
length = {r["name"]: int(r["length"]) if r["length"].strip().isdigit() else 0
          for r in csv.DictReader(open("panel/targets.csv"))}
rows = [l.split() for l in os.environ["TGT"].splitlines() if l.strip()]
load = [0] * n
lists = [[] for _ in range(n)]
for r in sorted(rows, key=lambda r: -length.get(r[0], 0)):
    i = min(range(n), key=lambda j: load[j])
    lists[i].append(" ".join(r))
    load[i] += length.get(r[0], 0)
for i, lst in enumerate(lists):
    with open(f"/tmp/panel_screen_w{i}.txt", "w") as fh:
        fh.write("\n".join(lst) + ("\n" if lst else ""))
    print(f"worker {i}: {[l.split()[0] for l in lst]} (load {load[i]} res)", file=sys.stderr)
EOF

    screen_worker() { # $1 worker list file
        while read -r name fs exp <&3; do
            y=$(yaml_of "$name")
            [ -f "$y" ] || { echo "[skip] $name: $y missing (run fetch)"; continue; }
            for s in $SEEDS; do
                if [ -s "eval_${name}_stock_s${s}.csv" ]; then
                    echo "[skip] $name seed $s: screen CSV exists"; continue
                fi
                echo "=== screen $name (seed $s, $N samples, gpu ${CUDA_VISIBLE_DEVICES:-default}) ==="
                boltz predict "$y" --use_msa_server --use_potentials < /dev/null \
                    --diffusion_samples "$N" --sampling_steps "$STEPS" $BOLTZ_EXTRA_ARGS \
                    --output_format pdb --seed "$s" --cache "$BOLTZ_CACHE" \
                    --out_dir "results_${name}_stock_s${s}" --override \
                    || { echo "[FAIL boltz] $name seed $s"; continue; }
                run_eval "results_${name}_stock_s${s}" "$name" "$fs" "stock_s${s}" "$(yaml_stem "$y")" \
                    || { echo "[FAIL eval] $name seed $s"; continue; }
                [ "$s" = "42" ] && freeze_msa "$name" "results_${name}_stock_s${s}"
            done
        done 3< "$1"
    }

    if [ "$N_GPUS" -le 1 ]; then
        screen_worker /tmp/panel_screen_w0.txt
    else
        PIDS=""
        for i in $(seq 0 $((N_GPUS-1))); do
            ( export CUDA_VISIBLE_DEVICES=$i; screen_worker "/tmp/panel_screen_w$i.txt" ) \
                > "screen_gpu$i.log" 2>&1 &
            PIDS="$PIDS $!"
            echo "[worker $i] GPU $i -> screen_gpu$i.log (pid $!)"
        done
        rc=0
        for p in $PIDS; do wait "$p" || rc=1; done
        for i in $(seq 0 $((N_GPUS-1))); do
            echo "[worker $i] failures: $(grep -c '^\[FAIL' "screen_gpu$i.log" 2>/dev/null || true) (screen_gpu$i.log)"
        done
        [ "$rc" -eq 0 ] || echo "WARNING: a worker exited non-zero; check screen_gpu*.log" >&2
    fi

    # Serial diagnosis pass, AFTER all workers: panel_diagnose rewrites
    # panel/panel_diagnosis.csv whole-file, so concurrent workers would
    # race; diagnosing here keeps the table deterministic.
    while read -r name fs exp <&3; do
        csvs=""; dirs=""; missing=0
        for s in $SEEDS; do
            csvs="$csvs eval_${name}_stock_s${s}.csv"
            dirs="$dirs $(pred_dir results_${name}_stock_s${s} "$name")"
            [ -s "eval_${name}_stock_s${s}.csv" ] || missing=1
        done
        if [ "$missing" -eq 1 ]; then
            echo "[incomplete] $name: seed CSV(s) missing - rerun screen to retry (log above / screen_gpu*.log)"
            continue
        fi
        python3 panel/panel_diagnose.py --name "$name" \
            --csvs $csvs --pred_dirs $dirs \
            --basin_rule "$(basin_rule "$name" "$fs")" \
            --expected "$([ "$exp" = "ref2" ] && echo closed || echo open)" \
            --out "diagnosis_${name}.txt"
    done 3<<< "$TARGETS"
    echo; echo "=== screen done; panel/panel_diagnosis.csv: ==="
    column -s, -t panel/panel_diagnosis.csv 2>/dev/null || cat panel/panel_diagnosis.csv
    ;;

intervene)
    patch_gate
    [ -f panel/panel_diagnosis.csv ] || { echo "run screen first" >&2; exit 1; }
    python3 - "$N_HEALTHY_CONTROLS" <<'EOF' > /tmp/panel_interventions.txt
import csv, random, sys
cap = int(sys.argv[1])
DEV = {"adk", "mbp"}
meta = {r["name"]: r for r in csv.DictReader(open("panel/targets.csv"))}
refrmsd = {}
try:
    for r in csv.DictReader(open("panel/fetch_report.csv")):
        try:
            refrmsd[r["name"]] = float(r["ref_rmsd"])
        except (ValueError, KeyError):
            pass
except FileNotFoundError:
    pass

def is_primary(n):
    rr = refrmsd.get(n)
    return (n not in DEV and meta.get(n, {}).get("set") == "oc23"
            and rr is not None and rr == rr and rr >= 4.0)

healthy_primary = []
for r in csv.DictReader(open("panel/panel_diagnosis.csv")):
    if r["verdict"] in ("damp", "boost") and r["alpha"]:
        contra = "2.5" if r["verdict"] == "damp" else "0.5"
        print(r["name"], "ctrl", r["alpha"])    # diagnosed arm
        print(r["name"], "contra", contra)      # contra-diagnosed arm (sign-law control)
    elif r["verdict"] == "healthy" and is_primary(r["name"]):
        healthy_primary.append(r["name"])

# Healthy controls (amendment 2026-10-05b, audit round 2 point 3): ALL
# primary-stratum healthy targets, up to cap; if more than cap, a
# fixed-seed (7) draw, logged. With <=5 controls this is a WEAK
# specificity claim and is worded as such in the paper.
pick = sorted(healthy_primary)
if len(pick) > cap:
    pick = sorted(random.Random(7).sample(pick, cap))
    note = f"fixed-seed(7) draw {pick} from {sorted(healthy_primary)}"
else:
    note = "all primary healthy targets" if pick else "none"
with open("panel/healthy_controls.txt", "w") as fh:
    fh.write(f"healthy controls ({note}), cap={cap}\n")
    for n in pick:
        fh.write(n + "\n")
print(f"# healthy controls ({note}): {pick}", file=sys.stderr)
for n in pick:
    print(n, "ctrl", "0.5")
    print(n, "contra", "2.5")
EOF
    while read -r name tag alpha <&3; do
        [ -n "$name" ] || continue
        if [ -n "$ONLY" ] && [[ ",$ONLY," != *",$name,"* ]]; then continue; fi
        fs=$(awk -F, -v n="$name" '$1==n {print $12}' panel/targets.csv)
        run_arm "$name" "$fs" "$tag" "$alpha"
    done 3< /tmp/panel_interventions.txt
    echo "intervene done. Next: ./panel/run_panel.sh analyze"
    ;;

# Pre-registered secondary arm (amendments 2026-10-05b/c): MSA-depth-
# engineered boost. Runs on ALL selected healthy controls (the same
# targets as the no-harm controls, from panel/healthy_controls.txt) —
# no per-target cherry-picking. Depth frozen at 32 in advance, never
# tuned per target. Manipulation check: the alpha=2.5 arm runs only if
# the DEGRADED stock (depth-32 MSA, no control) re-diagnoses as 'boost';
# pass rate reported by analyze.
degrade)
    patch_gate
    [ -f panel/panel_diagnosis.csv ] || { echo "run screen first" >&2; exit 1; }
    [ -s panel/healthy_controls.txt ] || { echo "no healthy controls; run intervene first" >&2; exit 1; }
    while read -r DNAME <&3; do
        [ -n "$DNAME" ] || continue
        [ -s "panel_msa/$DNAME.a3m" ] || { echo "[skip] $DNAME: no frozen MSA"; continue; }
        fs=$(awk -F, -v n="$DNAME" '$1==n {print $12}' panel/targets.csv)
        exp=$(awk -F, -v n="$DNAME" '$1==n {print $17}' panel/targets.csv)
        python3 panel/msa_subsample.py "panel_msa/$DNAME.a3m" \
            "panel_msa/$DNAME.d${DEGRADE_DEPTH}.a3m" --depth "$DEGRADE_DEPTH"
        sed "/      sequence:/a\\      msa: panel_msa/$DNAME.d${DEGRADE_DEPTH}.a3m" \
            "yamls/$DNAME.yaml" > "yamls/$DNAME.d${DEGRADE_DEPTH}.yaml"
        local_dstem="$DNAME.d${DEGRADE_DEPTH}"
        # 1) degraded STOCK (no control) + diagnosis = manipulation check
        csvs=""; dirs=""
        for s in $SEEDS; do
            csvs="$csvs eval_${DNAME}_degrade_stock_s${s}.csv"
            dirs="$dirs $(pred_dir results_${DNAME}_degrade_stock_s${s} "$local_dstem")"
            if [ -s "eval_${DNAME}_degrade_stock_s${s}.csv" ]; then
                echo "[skip] $DNAME degrade_stock seed $s: exists"; continue
            fi
            if [ -d "$(pred_dir results_${DNAME}_degrade_stock_s${s} "$local_dstem")" ]; then
                echo "[reuse] $DNAME degrade_stock seed $s: predictions exist, eval only"
            else
                echo "=== $DNAME degrade_stock (MSA depth $DEGRADE_DEPTH, no control) seed $s ==="
                boltz predict "yamls/$DNAME.d${DEGRADE_DEPTH}.yaml" --use_potentials < /dev/null \
                    --diffusion_samples "$N" --sampling_steps "$STEPS" $BOLTZ_EXTRA_ARGS \
                    --output_format pdb --seed "$s" --cache "$BOLTZ_CACHE" \
                    --out_dir "results_${DNAME}_degrade_stock_s${s}" --override \
                    || { echo "[FAIL boltz] $DNAME degrade_stock seed $s"; continue; }
            fi
            run_eval "results_${DNAME}_degrade_stock_s${s}" "$DNAME" "$fs" "degrade_stock_s${s}" "$local_dstem" \
                || echo "[FAIL eval] $DNAME degrade_stock seed $s"
        done
        dstock_ok=1
        for s in $SEEDS; do
            [ -s "eval_${DNAME}_degrade_stock_s${s}.csv" ] || dstock_ok=0
        done
        if [ "$dstock_ok" -eq 0 ]; then
            echo "[incomplete] $DNAME degrade_stock: CSV(s) missing - rerun degrade"
            continue
        fi
        python3 panel/panel_diagnose.py --name "$DNAME" \
            --csvs $csvs --pred_dirs $dirs \
            --basin_rule "$(basin_rule "$DNAME" "$fs")" \
            --expected "$([ "$exp" = "ref2" ] && echo closed || echo open)" \
            --out "diagnosis_${DNAME}_degrade.txt" \
            --report_csv panel/panel_diagnosis_degrade.csv
        dv=$(awk -F, -v n="$DNAME" '$1==n {print $8}' panel/panel_diagnosis_degrade.csv | tail -1)
        if [ "$dv" != "boost" ]; then
            echo "  $DNAME: degraded-stock diagnosis = $dv -> manipulation FAILED, no boost arm"
            continue
        fi
        echo "  $DNAME: degraded stock diagnoses as boost -> manipulation VALID"
        # 2) boost arm on the degraded background (same degraded MSA yaml)
        for s in $CTRL_SEEDS; do
            if [ -s "eval_${DNAME}_degrade_s${s}.csv" ]; then
                echo "[skip] $DNAME degrade seed $s: exists"; continue
            fi
            if [ -d "$(pred_dir results_${DNAME}_degrade_s${s} "$local_dstem")" ]; then
                echo "[reuse] $DNAME degrade seed $s: predictions exist, eval only"
            else
                echo "=== $DNAME degrade (MSA depth $DEGRADE_DEPTH) alpha=2.5 tw=0.5 seed $s ==="
                NU_ALPHA=2.5 NU_TW=0.5 \
                boltz predict "yamls/$DNAME.d${DEGRADE_DEPTH}.yaml" --use_potentials < /dev/null \
                    --diffusion_samples "$N" --sampling_steps "$STEPS" $BOLTZ_EXTRA_ARGS \
                    --output_format pdb --seed "$s" --cache "$BOLTZ_CACHE" \
                    --out_dir "results_${DNAME}_degrade_s${s}" --override \
                    || { echo "[FAIL boltz] $DNAME degrade seed $s"; continue; }
            fi
            run_eval "results_${DNAME}_degrade_s${s}" "$DNAME" "$fs" "degrade_s${s}" "$local_dstem" \
                || echo "[FAIL eval] $DNAME degrade seed $s"
        done
    done 3< <(grep -v '^#' panel/healthy_controls.txt | grep -v '^healthy controls')
    ;;

analyze)
    python3 panel/analyze_panel.py --targets panel/targets.csv --seeds $SEEDS --ctrl_seeds $CTRL_SEEDS
    ;;

*)
    echo "usage: $0 {fetch|screen|intervene|degrade|analyze} [--only a,b] [--set oc23|ms15|tp16]" >&2
    echo "  env knobs: GPUS=N (screen: split targets across N cards), BOLTZ_EXTRA_ARGS, N, SEEDS, CTRL_SEEDS" >&2
    exit 1 ;;
esac
