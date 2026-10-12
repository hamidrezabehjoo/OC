# panel/ — prospective panel screen (OC23 + MS15 [+ optional TP16])

Screen-first prospective test of the noise controller across the two-state
benchmarks of Suzuki & Amagasa (JCIM 2026). **Operations (setup, cluster
jobs, knobs, failure triage) are in the top-level `README.md`** — this file
documents the internals for method readers.

Everything is driven by `targets.csv`; the frozen pre-registration is
`prospective_test_panel_protocol.md` (v2.5, git tag `panel-v2.5`).

## Phase flow (all via `run_panel.sh`, all resumable)

| phase | what it does | key outputs |
|---|---|---|
| `fetch` | UniProt/RCSB download, chain extraction, Needleman–Wunsch renumbering onto the input sequence, yaml writer, sanity checks | `yamls/`, `panel_data/`, `fetch_report.csv` |
| `screen` | stock Boltz-2, 100 samples × 3 seeds (42, 1, 2) per target; first seed's MSA frozen to `panel_msa/` and reused by every later arm of that target; frozen diagnosis per target | `eval_*_stock_s*.csv`, `diagnosis_*.txt`, `panel_diagnosis.csv` |
| `intervene` | diagnosed damp/boost targets get the diagnosed arm AND the contra-diagnosed arm (sign-law control), 5 seeds each; all primary-stratum healthy targets (cap 5) get both as no-harm controls | `eval_*_{ctrl,contra}_s*.csv`, `healthy_controls.txt` |
| `degrade` | engineered boost on the same healthy controls: MSA subsampled to depth 32, degraded stock must re-diagnose as boost (manipulation check), then α = 2.5 arm | `panel_msa/*.d32.a3m`, `panel_diagnosis_degrade.csv` |
| `analyze` | per-target success (direction in all arm seeds + seed-level bootstrap CI), paired sign-law claim, efficacy claim, D1 geometry/reachability diagnostics | `panel_summary.csv` |

## Files

| file | role |
|---|---|
| `targets.csv` | target registry: refs, UniProt, labels, flags. One row per target. |
| `fetch_panel.py` | reference download + renumbering; NMR entries use the medoid model |
| `eval_panel.py` | residue-number intersection matching (references with missing density OK); Kabsch/TM/basin rules identical to the validated single-target pipeline. `--basin_rule tm` for fold-switchers; `--escape_rmsd auto` = min(3.5, 0.5 × ref-pair RMSD) |
| `panel_diagnose.py` | frozen decision table incl. ceiling clause and borderline guard (±3.0 pts) |
| `analyze_panel.py` | exact Wilcoxon, seed-level bootstrap, Fisher pooled (descriptive), Holm across basin rules if D1 fires |
| `msa_subsample.py` | deterministic depth-32 subsample of a frozen .a3m (engineered boost) |
| `make_manifest.py` | write/verify `MANIFEST.md5` — on the GPU machine use `--verify`, never re-fetch |
| `run_panel.sh` | the driver |
| `REGISTRATION.md` | tag/commit/DOI record (unpinned on purpose) |
| `outcome_paragraphs.md` | pre-written results paragraphs (Cases A/B/C + limitations) |
| `tm_check_*.csv` | TM-score reconciliation vs the AFsample2 Zenodo references (protocol §7) |

## Filling the remaining MS15 rows

Four rows in `targets.csv` are `status=tbd` pending Suzuki Table 2 values
(JCIM 2026, doi:10.1021/acs.jcim.6c02094). Fill
`ref1_pdb/ref1_chain/ref1_label/ref2_pdb/ref2_chain/ref2_label`, set
`status=ready`, then `fetch` + `screen` as usual. This is a dated protocol
amendment when it happens.

TP16 (16 transporters, optional phase 2):
`python3 panel/fetch_panel.py --targets panel/targets.csv --include-optional`,
then `./panel/run_panel.sh screen --set tp16`.
