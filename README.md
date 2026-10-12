# DiffusionControl — prospective panel test of inference-time noise control

Companion code for *Optimal Inference-Time Noise Control for Diffusion
Models*. The pipeline screens a frozen panel of two-state proteins with
**stock Boltz-2**, diagnoses each ensemble with a frozen decision table,
then runs the **frozen intervention arms** (damp α = 0.5 / boost α = 2.5,
window tw = 0.5) plus contra-diagnosed controls. The scientific design is
frozen in `panel/prospective_test_panel_protocol.md` (v2.5); this README
is the operator's manual.

## The four commands

Everything runs through one driver, from the repo root. Every phase is
**resumable**: finished work is skipped, only missing pieces run.

```bash
./panel/run_panel.sh fetch        # once: sequences, references, yamls (DONE — pinned)
./panel/run_panel.sh screen       # stock: 100 samples × 3 seeds per target + diagnosis
./panel/run_panel.sh intervene    # ctrl + contra arms on diagnosed targets, 5 seeds each
./panel/run_panel.sh degrade      # engineered-boost arm on healthy controls
./panel/run_panel.sh analyze      # panel-level statistics → panel/panel_summary.csv
```

Useful filters for any phase: `--only adk,rfah` (named targets),
`--set oc23|ms15` (stratum).

## Repo map

```
diffusionv2.py     the controller: stock boltz 2.2.1 sampler + NU_ALPHA hook
                   + SVD numerics guard. PINNED. Installed over boltz's copy
                   at setup (setup_env.sh does this).
setup_env.sh       one-time environment build (boltz[cuda]==2.2.1, patch,
                   verification gates)
panel/
  run_panel.sh     the driver — the only script you run
  targets.csv      the frozen 25-target registry
  *.py             fetch / eval / diagnose / analyze / manifest / MSA subsample
  prospective_test_panel_protocol.md   the frozen pre-registration (v2.5)
  MANIFEST.md5     md5s of every pinned file — the freeze gate
  REGISTRATION.md  tag / commit / DOI record (deliberately unpinned)
  outcome_paragraphs.md              pre-written results text (Cases A/B/C)
  README.md        panel internals for method readers
panel_data/        reference PDBs (pinned)
yamls/             boltz inputs (pinned)
```

Generated at runtime (never committed, never pinned — safe to keep):
`results_*`, `eval_*.csv`, `diagnosis_*.txt`, `screen_gpu*.log`,
`panel_msa/`, `yamls/*.msa.yaml`, `yamls/*.d32.yaml`,
`panel/panel_diagnosis.csv`, `panel/healthy_controls.txt`.

## Setup (once per machine)

```bash
bash setup_env.sh                 # NSCC ASPIRE 2A runbook: scratch env,
                                  # cu128 torch before boltz, nvrtc fix,
                                  # controller patch — idempotent, re-run
                                  # after any boltz reinstall/upgrade
conda activate /scratch/users/industry/ism/behjoo/envs/boltz_ctrl
python3 panel/make_manifest.py --verify   # must print: manifest OK
```

Scratch is purged after 30 days — archive `eval_*.csv`, `panel/panel_diagnosis.csv`,
`panel_msa/` and `panel/panel_summary.csv` to durable storage after every phase.

If you edit nothing, `--verify` always passes. **After every file transfer
to the GPU machine, run `--verify` first** — it catches truncated or stale
copies instantly. Never re-run `fetch` on the GPU machine (it would
re-download references and break the freeze).

## Running on the cluster (PBS)

`submit_panel.sh` (repo root) is the batch job. Submit from `~/DiffusionControl`:

```bash
qsub -v PHASE=screen     submit_panel.sh
qsub -v PHASE=intervene  submit_panel.sh    # your current phase
qsub -v PHASE=degrade    submit_panel.sh
qsub -v PHASE=analyze    submit_panel.sh
qsub -v PHASE=all        submit_panel.sh    # everything, in order
```

It activates the scratch env, sets the nvrtc/caches exports (batch jobs don't
read `~/.bashrc`), auto-detects `GPUS` from the allocation, defaults
`BOLTZ_EXTRA_ARGS="--max_parallel_samples 30"`, runs `make_manifest.py
--verify` before any compute, and chains phases for `all`. Resources are in
the `#PBS` header (`select=1:ncpus=8:ngpus=2`, 24 h) — edit there if your
allocation differs. Every phase is resumable (`[skip]` on finished seeds), so
a walltime kill means: resubmit the same PHASE.

## Environment knobs

| knob | default | meaning |
|---|---|---|
| `GPUS` | 1 | screen only: split targets across N GPUs (one boltz per card, length-balanced). boltz never pools GPU memory — this buys time, not memory. |
| `BOLTZ_EXTRA_ARGS` | empty | extra flags for every `boltz predict`. For proteins ≳450 aa on 40 GB cards use `--max_parallel_samples 30` (≤10 samples/forward) or `34` (≤4). **Trap:** with 100 samples, values 5/10/20/25/50 are silently no-ops (boltz chunk formula: `100 % m + 1`). Provably does not change per-sample noise streams — execution-only. |
| `N`, `STEPS`, `SEEDS`, `CTRL_SEEDS` | 100, 200, "42 1 2", "42 1 2 7 13" | **frozen by the protocol. Do not change.** |

## Log lines and what to do

| line | meaning | action |
|---|---|---|
| `[skip] <name> seed s: ...` | already done | none — resume working as intended |
| `[reuse] <name> <tag> seed s: predictions exist` | structures on disk from a previous run; only the deterministic eval re-runs | none — this is how the 2026-10-08 eval-path fix recovers completed runs |
| `CORRUPT: panel_msa/<name>.a3m has NUL bytes` | frozen MSA damaged (boltz's writer appends a trailing NUL; fixed at freeze time for new freezes) | strip the byte: `python3 -c "p='panel_msa/<name>.a3m'; d=open(p,'rb').read(); assert d[-1:]==b'\x00' and d.count(b'\x00')==1; open(p,'wb').write(d[:-1])"`, log the new md5, rerun |
| `[worker i] GPU i -> screen_gpu<i>.log` | screen split across GPUs | progress is in those logs |
| `WARNING: ran out of memory` | GPU OOM (protein too big for unchunked sampling) | set `BOLTZ_EXTRA_ARGS="--max_parallel_samples 30"` and rerun; if still OOM, `34` |
| `[nu-guard] gesvd failed...` | SVD fallback engaged (rare, boost arms) | none — counted for the paper; many per run = tell us |
| `tarfile.ReadError` (MSA server) | colabfold returned a corrupt tarball | transient — rerun |
| `[incomplete] <name>` | some seeds failed this pass | fix cause above, rerun same command |
| `MSA_FLAG: unbound variable` / silent exit | stale or truncated `run_panel.sh` | redeploy, `--verify` |

## The freeze (what you may and may not touch)

Pinned by `panel/MANIFEST.md5` — never edit on the GPU machine:
`diffusionv2.py`, `panel/*.py`, `panel/run_panel.sh`, `targets.csv`, the
protocol, `panel_data/`, `yamls/`, and `panel_msa/*.a3m` once frozen.

Change flow for pinned files: edit in the source repo →
`make_manifest.py` → commit → move tag `panel-v2.5` → deploy →
`--verify` on the cluster. External registration (GitHub release →
Zenodo DOI) is recorded in `panel/REGISTRATION.md`, which is unpinned on
purpose.
