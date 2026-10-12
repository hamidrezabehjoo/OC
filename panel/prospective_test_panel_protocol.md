# Prospective panel screen — frozen protocol

**Version 2.5, 2026-10-05 (e).** Panel-level prospective test (advisor
comment 10 / screen-first design). All constants identical to v1 (single-
target MBP protocol) unless an amendment below says otherwise. Frozen as
git tag `panel-v2.5` with `panel/MANIFEST.md5` pinning every input file.
**External timestamp (required before the screen run)**: push tag
`panel-v2.5` to the project repo and mint a Zenodo DOI via a GitHub
release (or an OSF registration); a local tag is not a timestamp. The
DOI itself is recorded in `panel/REGISTRATION.md`, which is deliberately
**not** pinned by the manifest (pinning it would change the md5 after
archival and invalidate the tag it documents). When moving this tree
into the real project repo: commit the same content, run
`python3 panel/make_manifest.py --verify`, and only then tag and
release — the Zenodo record must come from that repo.

**Benchmark version.** We use the JCIM 2026 organization (OC23 + MS15
[+ TP16]; Suzuki & Amagasa, doi:10.1021/acs.jcim.6c02094). The bioRxiv v2
text describes a differently assembled 39-soluble + 47-transporter set
(BioEmu + OC23 merged; IOMemP + Wu & Feng). The JCIM table is the cited
benchmark; OC23 rows come verbatim from AFsample2 Table S1
(doi:10.1038/s42003-025-07791-9).

## 0. Development cases (NOT in the primary tally)

| id | target | role | outcome |
|----|--------|------|---------|
| A0 | ADK | the diagnosis rule, gain α=0.5, and the success endpoint were calibrated/validated here | published damp success: open 15%→29%, p=4.8e-05 (pooled 3 seeds × 100) |

**ADK anchor note (2026-10-05(c)).** Round-1 drafts quoted the anchor as
18%→29%, p=0.035 — that was **seed 42 alone** (18/100 open), from the
single-seed screen design. The 2026-10-05 three-seed amendment pools
45/165/90 = 15.0% stock open (the paper's published number; an
independent rerun reproduced 15.3%), hence p tightens to 4.8e-05. The
basin definition never changed: the paper's `eval_basins.py` uses
escape_rmsd = 3.5 Å (default) and `eval_panel.py --escape_rmsd auto`
gives min(3.5, 0.5 × 7.13) = min(3.5, 3.565) = **3.5 Å** — identical
threshold and identical assignment rule (`r ≤ esc and r < r_other`),
verified line-by-line. The radius cap is not the cause of the 18→15
move; seed pooling is.
| M0 | MBP | ceiling clause written after its stock diagnosis | ceiling (open 100%, pairwise 0.30 Å), no intervention |

Both are development cases: they shaped the protocol, so they cannot be
evidence for it. Reported separately, always labeled "development".
The ADK anchor check (diagnose → damp, α=0.5 → open 15%→29%, p<0.05) is
a **pipeline regression check** — the protocol was designed to pass it —
not validation of the method (amendment 2026-10-05(b), audit round 2,
point 6); it is never cited as evidence.

## 1. Panel and strata

- **oc23-large** (PRIMARY): OC23 targets with reference-pair Cα RMSD
  ≥ 4.0 Å. Count: **12** (recomputed 2026-10-05(b) after the NMR-medoid
  amendment; boundary-fragile rows: q18a65 4.03, p00558 3.91, a0qtt2 4.14).
  Membership frozen by name (2026-10-05(c)) — the stratum cannot drift:
  a0a075q0w3, a0qtt2, a2rj53, a6uvt1, o76728, p21589, p31133, p40131,
  p62495, q18a65, q5f9m1, q9ere7.
- **oc23-xl** (robustness subset of oc23-large; 2026-10-05(c)):
  ref-pair RMSD ≥ 5.0 Å — a2rj53 (5.07), o76728 (5.02), p21589 (11.08),
  p40131 (10.99), q5f9m1 (5.75), q9ere7 (20.62), a0a075q0w3 (7.47).
  Reported as a free robustness check on the primary endpoint (the 4.0 Å
  boundary sits only 0.23 Å above the MBP calibration point, and 5 of the
  12 primaries lie in 4.0–5.0 Å).
- **oc23-small** (secondary): OC23 targets below 4.0 Å — basin radius
  approaches model error; analyzed with the scaled threshold (§3) and
  reported as a separate stratum. Count: **11**.
- The 4.0 Å cutoff is measured on the **same metric** as the MBP
  development case (3.77 Å): all-Cα Kabsch RMSD on shared positions,
  same code path (`fetch_panel.py`), so the stratum boundary is
  internally consistent with the calibration point.
- **ms15** (exploratory): hand-picked rows (anchor ADK, RfaH; remaining
  rows marked `tbd` until filled from JCIM Table 2 — a dated amendment,
  before those rows are screened).
- **tp16** (optional phase 2): 16 transporters, `optional=1`; SPF1
  excluded. MurJ (b7ie18) also appears in OC23; tm_0287/tm_0287v2 share
  ref 6QV1 — deduplicate if TP16 is ever pooled with OC23.
- Strata are never pooled: a fold-switcher failure must not read as a
  failure of the hinge-motion claim, and vice versa.
- **Expected state** (frozen 2026-10-05, wording fixed 2026-10-05(d)):
  `targets.csv` column 17 (`expected_state`) is **ref1 for every row**,
  where ref1 is the apo/ground state. For RfaH that means the α-hairpin
  (2OUG) is expected and the β-barrel (6C6S) is the activated minority
  state. The diagnosis is therefore **reference-informed under a
  two-state hypothesis**, not reference-free: basins are defined by
  reference locations and the boost branch needs to know which state is
  the wrong one. This is how the contrast with Suzuki et al.'s
  target-specific β sign should be read — they tune a per-target scalar
  informed by known states; we act on a per-target *diagnosis* informed
  by the same two-state references, with the sign of the intervention
  fixed by the diagnosis rather than tuned.

## 2. Frozen decision table (per target, on the pooled stock screen)

a. intermediate fraction ≥ 0.50 → **damp** (α = 0.5)
b. else one basin share ≥ 0.50 AND mean pairwise Cα RMSD < 2.0 Å →
   collapsed: dominant = expected state with share ≥ 0.95 → **ceiling**
   (no run; amendment 2026-10-01); otherwise → **boost** (α = 2.5)
c. otherwise → **healthy** (no intervention)
d. **indeterminate** (amendment 2026-10-05; guard widened 2026-10-05(b)):
   point estimate within **3.0 points** of a decision threshold — ≥ 1
   binomial SE at n = 300, p = 0.5 (SE = 2.9 pt) — OR per-seed
   point-estimate verdicts disagree → no intervention; reported as its own
   bin. Expected noise-driven verdict flips at n = 300: a true fraction
   2.5 pt from the threshold crosses it ~19% of the time; 5 pt out, ~4%.
   Bootstrap CIs of the decision statistics are reported for transparency
   but do not gate (at practical n they straddle 0.50 for exactly the
   near-threshold targets the panel cares about; the guard targets
   coin-flip decisions).

## 3. Screen settings (frozen)

- Boltz-2 v2.2.1, `--use_potentials`, 200 sampling steps, 100 diffusion
  samples × **3 seeds (42, 1, 2)** per target (amended from 1×100 on
  2026-10-05 for diagnosis stability; matches the ADK/MBP design).
- Input sequence: construct-level, ref1 SEQRES, uniformly for all rows
  (amendment 2026-10-05; one input policy, no UniProt/construct mixing).
- **MSA freezing**: the first stock seed runs with `--use_msa_server`;
  the resulting .a3m is copied to `panel_msa/<name>.a3m` (md5 logged) and
  reused by every later arm of that target. No cross-arm MSA drift.
- Basins: RMSD rule, basin radius = **min(3.5, 0.5 × ref-pair Cα RMSD)**
  (amendment 2026-10-04(b)). The 3.5 Å cap is hard: even for the widest
  hinge in the panel (q9ere7, ref-pair 20.6 Å) the basin radius is
  3.5 Å, not 10 Å — no target gets an unbounded basin (audit round 2,
  point 1). Fold-switch rows (`foldswitch=1`): TM rule (TM ≥ 0.5, margin
  ≥ 0.15).
- **Sensitivity analysis** (pre-registered 2026-10-05(b), d0 fixed
  2026-10-05(e)): every diagnosis and every success call is additionally
  recomputed under the TM basin rule (`--basin_rule tm`) with the
  **standard length-dependent d0** (1.24·(L−15)^⅓ − 1.8) — the fixed
  d0 = 3.5 Å convention exists only to match AFsample2's table (their
  TMalign runs at `-d 3.5`) and is far too strict for 700-residue
  targets; it is never used for the sensitivity re-diagnosis.
  Discrepancies between the two rules are reported per target. The RMSD
  rule remains the primary one. **Multiplicity**: if D1 fires and the
  TM-rule analysis becomes co-primary, the sign-law claim is
  **Holm-corrected across the two basin-rule analyses** (two tests;
  2026-10-05(e)); analyze_panel computes both p-values and the Holm
  verdict unconditionally.
- **NMR references** (2026-10-05(b), pre-registered): for multi-model
  (NMR) entries the **medoid model** — the model minimizing mean pairwise
  Cα RMSD to all other models, ties broken by lowest model number — is
  the reference, replacing the earlier MODEL-1 rule. Affects q9ere7
  (2KX2-class entries) and p62495; both stay in oc23-large.
- **Expected state**: frozen per row in `targets.csv` (`expected_state`
  column, dated 2026-10-05; ref1 = apo/expected for all current rows).
  The diagnosis uses the stock ensemble plus this frozen label; it never
  uses controlled-arm data.

## 4. Interventions (frozen arms)

- Diagnosed targets (damp/boost): the **diagnosed arm** (α=0.5 or 2.5,
  tw=0.5) AND the **contra-diagnosed arm** (the other α) — the sign-law
  control (amendment 2026-10-05). Prediction: the contra arm does not
  improve the endpoint; if it does, the benefit is generic diversification,
  not the sign law.
- **Healthy controls** (amended 2026-10-05(b)): ALL healthy-diagnosed
  targets **in the primary stratum**, up to 5, receive both arms; if more
  than 5 qualify, a fixed-seed (7) random draw, logged with the draw list
  to `panel/healthy_controls.txt`. (The earlier "first 3 in panel order"
  rule was biased: registry order correlates with motion size.) With
  ≤ 5 controls this is a **weak specificity claim** — worded as such in
  the paper. Criterion: no significant *decrease* of expected-state
  occupancy under either arm ("the method knows when not to act" is
  tested, not asserted).
- **Engineered boost arm** (2026-10-05(b), revised 2026-10-05(c),
  secondary): the boost bin is expected to have < 4 members and thus be
  inconclusive; to give the boost regime guaranteed data points, the
  **healthy-control targets** (the same ≤ 5 primary-stratum targets as
  the no-harm controls — no separate selection, no order-dependence) are
  each run with their frozen MSA deterministically subsampled to
  **depth 32, fixed in advance and never tuned per target**
  (`panel/msa_subsample.py`: query + constant stride, no RNG).
  **Manipulation check**: first a degraded *stock* run (depth-32 MSA, no
  control) is diagnosed; the α = 2.5 arm proceeds only where this
  degraded stock re-diagnoses as **boost** (collapsed into one basin),
  and the pass rate is reported. Scoring is against the degraded stock,
  by the boost success criterion; reported as a secondary arm, never
  pooled with the spontaneous boost bin.
- ceiling / indeterminate: no run.
- **5 seeds × 100 per arm** (42, 1, 2, 7, 13; amended 2026-10-05(e)):
  at 3 seeds the seed-level bootstrap is degenerate — the 2.5th
  percentile equals the smallest seed shift (probability (1/3)³ = 3.7%
  > 2.5%), so "CI excludes 0" is exactly "3/3 positive" and the null
  rate is 1/2³ = 12.5%. At 5 seeds the criterion (5/5 direction +
  bootstrap CI) has null 1/2⁵ = 3.1%. If compute ever forces 3-seed
  arms, the paper must state plainly that the criterion is 3/3 with a
  12.5% per-target null rate. The stock screen stays at 3 seeds.
- Patch gate (NU_ALPHA hook) before any run.

## 5. Per-target success criteria (frozen)

- **damp**: expected-state occupancy ↑ vs stock, direction in ALL arm
  seeds (5/5) AND seed-level block-bootstrap 95% CI of the shift
  excludes 0 (amendments 2026-10-05(d)/(e); null rate 3.1% at 5 seeds).
  (ADK-calibrated endpoint: open 15%→29%.)
- **boost**: minority-basin occupancy ↑, same criteria.
- **Pseudo-replication rule** (2026-10-05(d)): predictions within a seed
  share the trunk and MSA, so pooled models are NOT independent.
  Per-seed shifts are reported for every arm; inference uses per-seed
  consistency + the seed-level bootstrap CI (resample seeds, 2000
  replicates). The pooled Fisher p is **descriptive only**.
- **healthy control**: expected-state occupancy shows no significant
  decrease (one-sided p ≥ 0.05) under either arm.

## 6. Panel-level endpoints

- **Two panel claims, fixed roles** (2026-10-05(c), roles fixed
  2026-10-05(d)):
  - **Efficacy claim**: pooled over ALL diagnosed **oc23-large** targets
    (damp + boost together), per-target ctrl-arm shifts signed in the
    diagnosed direction, **one-sided** exact Wilcoxon signed-rank test
    (two-sided cannot reach 0.05 below n = 6 — min p is 0.125 at n = 4,
    0.0625 at n = 5 — one-sided reaches 0.03125 at n = 5). This test can
    pass from generic broadening alone, so it is NOT the sign-law claim.
  - **Sign-law claim** (the primary scientific claim): **paired**
    per-target ctrl-vs-contra one-sided exact Wilcoxon — the diagnosed
    arm must beat its own contra arm. Reported with the contra
    false-positive rate (Wilson CI).
  - Per-target success fraction (Wilson 95% CI) reported alongside;
    damp / boost branches reported **descriptively** (success counts and
    median shifts — no per-branch binomial claims at these bin sizes).
- **Secondary**: per-stratum success (oc23-small, ms15, tp16), diagnosis-
  bin distribution, median occupancy shift, healthy-control no-harm rate,
  the **oc23-xl** robustness subset (§1), engineered-boost pass rate.
- **Pre-registered diagnostic D1** (2026-10-05(c), extended
  2026-10-05(d)): because the basin radius is capped at 3.5 Å,
  radius/separation is pinned at exactly 0.5 for every target with
  separation ≤ 7 Å (~19 of 25) and varies only for the rest, so the
  Spearman on radius/sep is underpowered (≈ 6 varying points). D1
  therefore has three parts:
  1. Spearman(intermediate fraction, radius/sep) and (…, length) —
     reported but underpowered by construction;
  2. **Group comparison**: median intermediate fraction for
     separation > 7 Å vs ≤ 7 Å, exact two-sided permutation test on
     medians — flagged a confound if p < 0.05 AND the median gap is
     > 0.10;
  3. **Reachability check**: best-of-300 stock RMSD to each reference vs
     the basin radius. A state the stock ensemble never reaches cannot
     host a diagnosable basin. Targets with BOTH references unreachable
     are annotated and any damp/boost diagnosis on them is treated as
     geometric (reported with the indeterminate bin); targets with one
     unreachable reference are annotated per row.
  **Conclusion rule**: if the group comparison flags, declare a geometry
  confound in the paper, report the diagnosis × geometry interaction,
  and the pre-registered TM-rule sensitivity diagnoses (§3) become
  co-primary for diagnosis; if also |ρ(length)| < 0.5, attribute the
  confound to basin geometry rather than size.
- **Power rule** (amended 2026-10-05(d)): a stratum/bin with < **5**
  members is reported as *inconclusive*, not failed (one-sided exact
  Wilcoxon needs n ≥ 5 to reach 0.05). In particular, if the boost bin
  ends up with 2–4 members, the boost branch is inconclusive — no
  failure claim.
- Literature context (not head-to-head — our per-target budget is 300
  stock models vs their 250-model/AF3-MSA protocol): vanilla Boltz-2
  dual-state coverage 43% (OC23), 33% (MS15), 25% (TP16).

## 7. Validation gates and exclusions

- **Reference gate** (per row, at fetch; amendment 2026-10-05, revised
  2026-10-05(b)/(c)/(d)): alignment identity ≥ 0.95 to the input sequence
  (hard gate — catches wrong chain/offset/transcription errors; NMR
  entries use the medoid model, §3). NMR caveat for the paper
  (2026-10-05(e)): q9ere7 and p62495 use OUR medoid model, not
  AFsample2's model choice, so comparisons to AFsample2/Suzuki numbers
  on those two targets are not like-for-like.
  **Settled against the source, 2026-10-05(d)**
  (`panel/tm_check_zenodo_2026-10-05.csv`): we downloaded AFsample2's
  actual reference files (Zenodo 10.5281/zenodo.14534088,
  `input_datasets.tar.zst` → `filtered_dict.pickle` + `pdbs/`). All 23
  OC23 reference pairs in our `targets.csv` are **exactly** the
  PDB:chain pairs AFsample2 used — our panel matches the OC23 definition.
  (The round-2/3 hypothesis that mismatched rows used different
  structures was wrong and is retracted.) The table's TM/RMSD columns
  come from the **TMalign binary run with a fixed d0 = 3.5 Å** and
  structure-based alignment (their `src/analyse_models.py`,
  `tmalign_and_extract_scores`), not from resseq-matched Kabsch: with
  d0 = 3.5 our metrics reproduce their TM within 0.06 for 37/45 targets;
  the 8 outliers are artifacts of naive resseq intersection (offset
  numbering, e.g. P15291 2fy7_A 122-398 vs 2fyc_B 131-402; NMR model
  choice for Q9ERE7), not different structures. Our pipeline is immune
  because it renumbers by sequence alignment. Table values remain
  informational; the identity gate stands.
- **Freeze verification** (2026-10-05(c)): the panel is frozen as git
  tag `panel-v2.3`; `panel/MANIFEST.md5` pins md5s of targets.csv, panel
  scripts, this protocol, prepared references, yamls and frozen MSAs.
  On the GPU machine: run `python3 panel/make_manifest.py --verify` —
  do NOT re-run fetch there (a re-fetch hits RCSB again and could change
  references).
- fetch failure; Boltz failure after one retry; < 80 usable screen models.
- `fetch_report.csv` and `panel_diagnosis.csv` keep the **latest row per
  target** (rewritten, not appended).

## 8. Cost estimate (A100-class)

Screen: 23 OC23 + 2 ms15 = 25 targets × 300 ≈ 7,500 predictions
≈ 13 GPU-hours. Interventions: diagnosed targets × 600 (both arms) +
≤ 5 healthy × 600 + 1 degrade arm × 300. TP16 optional phase 2:
15 targets, ~3–5× cost.

## Amendments

- **2026-10-01** (from v1): ceiling clause (§2) after the MBP stock
  diagnosis (open = 100%); recorded before any controlled MBP run.
- **2026-10-04**: single-target → panel-level. damp endpoint fixed to
  expected-state occupancy (ADK anchor: significant at p = 0.035; the
  intermediate-decrease variant fails on the anchor at p = 0.30 and was
  rejected before any panel run).
- **2026-10-04 (b)**: degenerate reference pairs (< 3.5 Å) not excluded —
  basin threshold auto-scales to min(3.5, 0.5 × ref RMSD); low-coverage
  rows rescued at construct level. (Superseded in part on 2026-10-05 by
  the oc23-large/oc23-small split.)
- **2026-10-05** (external audit round): construct-level input made
  uniform; NMR MODEL-1 parsing fix (2RQK-class entries); targets.csv
  parse fixes (eftu field count, spf1 quoting); `expected_state` column
  frozen; reference identity gate + TM review flag (RMSD-vs-table gate
  rejected as convention-incomparable); MSA freezing; 3-seed screen;
  borderline/per-seed-agreement indeterminate bin; contra-diagnosed arms
  and healthy controls added; strata defined (oc23-large primary);
  ADK/MBP reclassified as development cases; power/stopping rule;
  JCIM version cited as the benchmark definition.
- **2026-10-05 (b)** (external audit, round 2): (1) basin-radius cap
  stated explicitly (radius = min(3.5, 0.5 × ref RMSD) — the
  ~10 Å-basin scenario cannot occur) and TM-rule sensitivity analysis
  pre-registered for all targets; confirmed the 4.0 Å stratum cutoff and
  MBP's 3.77 Å share one metric/code path; (2) NMR references use the
  medoid model (pre-registered), refs rebuilt, strata recomputed —
  12/11 split confirmed, boundary rows flagged; (3) healthy controls =
  all primary-stratum healthy targets up to 5, else fixed-seed(7) draw
  logged; claim worded as weak specificity; (4) primary sign-law
  comparison = paired per-target ctrl-vs-contra (exact Wilcoxon) +
  contra false-positive rate with Wilson CI; (5) borderline guard widened
  2.5 → 3.0 pt (≥ 1 SE at n = 300), expected flip rates stated in §2;
  (6) ADK anchor relabeled a pipeline regression check, not validation;
  (7) MSA-depth-engineered boost arm (depth 32, deterministic stride
  subsample of the frozen MSA, one logged primary-healthy target)
  pre-registered as secondary; (8) table RMSD/TM verified on all 25
  targets and downgraded to informational (non-systematic deviations;
  identity gate 1.00 everywhere stands as the hard guard).
- **2026-10-05 (d)** (external audit, round 4): (1) expected-state freeze
  confirmed (column 17, all rows ref1) and semantics documented incl.
  RfaH; claims reworded "reference-informed under a two-state
  hypothesis"; (2) D1 extended: group comparison (sep ≤ 7 vs > 7 Å,
  exact permutation test) + reachability check (best-of-300 RMSD vs
  radius; both-unreachable ⇒ geometric diagnosis, reported with
  indeterminate); (3) test roles fixed: paired ctrl-vs-contra = the
  sign-law claim, pooled-vs-zero = efficacy claim; MIN_BIN 4 → 5;
  (4) TM question settled against AFsample2's Zenodo files: all 23
  reference pairs match exactly; their TM column is TMalign with fixed
  d0 = 3.5 Å (37/45 reproduced within 0.06; outliers are naive
  resseq-intersection artifacts); round-2/3 "different structures"
  hypothesis retracted; (5) pseudo-replication fix: per-seed shifts +
  seed-level block-bootstrap CI gate success; pooled Fisher p demoted to
  descriptive; (6) external timestamp made a gate for the screen run:
  push tag + GitHub release → Zenodo DOI (or OSF) before screening.
- **2026-10-05 (e)** (external audit, round 5): (1) DOI record moved to
  unpinned `panel/REGISTRATION.md` (tag, commit hash, md5 of the
  manifest) — pinning it would invalidate the tag it documents;
  (2) controlled arms go to **5 seeds** (null 3.1%; the 3-seed bootstrap
  is degenerate — "CI excludes 0" ≡ "3/3 positive", null 12.5%);
  (3) TM-rule sensitivity uses standard length-dependent d0 — fixed
  d0 = 3.5 Å was only for matching AFsample2's table; (4) Holm
  correction pre-registered across the two basin-rule analyses if D1
  fires; NMR medoid caveat (q9ere7, p62495) noted for the paper;
  (5) outcome paragraphs pre-written (`panel/outcome_paragraphs.md`),
  including the limitation that expected_state = ref1 means success =
  recovery of the deposited apo state, not true solution populations.
- **2026-10-05 (c)** (external audit, round 3): (1) ADK anchor move
  explained and documented (§0): 18% was seed-42-only under the
  single-seed design; pooled 3-seed 15.0% is the paper's number;
  basin definitions verified identical (3.5 Å, same rule) — radius cap
  innocent; (2) geometry diagnostic D1 pre-registered with a conclusion
  rule (§6); (3) primary test redefined: pooled signed shifts, one-sided
  exact Wilcoxon (two-sided can't reach 0.05 below n = 6); damp/boost
  branches descriptive; paired contrast one-sided; (4) ΔTM listed for
  all 23 OC23 rows + chain-level scan on the worst six — no chain
  combination reproduces the table, references are ours alone
  (`tm_check_2026-10-05.csv`); (5) engineered boost runs on ALL healthy
  controls, depth 32 frozen, manipulation check = degraded stock must
  re-diagnose as boost, pass rate reported; (6) git tag `panel-v2.3` +
  md5 manifest (verify on the GPU machine, never re-fetch), primary
  stratum membership frozen by name, oc23-xl ≥ 5 Å robustness subset
  pre-registered.
- **2026-10-06** (execution-only, no protocol change): `run_panel.sh`
  gains (a) a `BOLTZ_EXTRA_ARGS` env knob (default empty) so GPU-memory
  batching (`boltz predict --max_parallel_samples m`) can be set from
  the job scheduler without editing a pinned file, and (b) a `GPUS=N`
  knob that partitions screen targets across N cards (one boltz process
  per GPU, greedy length-balanced assignment), with diagnosis moved to
  a serial pass after all workers finish (panel_diagnose rewrites
  panel_diagnosis.csv whole-file, so concurrent diagnosis would race).
  boltz never pools memory across GPUs, so (b) does not remove the need
  for (a). Verified against boltz 2.0.0/2.1.1/2.2.1 source
  (`model/modules/diffusionv2.py`): initial noise, per-step noise, and
  random augmentation are drawn for the full multiplicity OUTSIDE the
  sample-chunk loop; chunking splits only the deterministic denoiser
  forward, so per-sample noise streams and all frozen settings
  (100 samples/seed, seeds, 200 steps) are unchanged. Note: with 100
  samples, m ∈ {5 (default), 10, 20, 25, 50} is a silent no-op
  (chunk count = 100 % m + 1 = 1); valid values are m = 30
  (≤ 10 samples/forward) or m = 34 (≤ 4 samples/forward).
- **2026-10-07** (execution-only, no protocol change): (1) The first
  boost (α = 2.5) arms exposed a numerical failure: cuSOLVER's fp32
  gesvd in boltz's per-step `weighted_rigid_align` can fail to converge
  on near-singular covariance matrices (denoiser output momentarily
  degenerate under boosted noise), killing the whole 100-sample run
  deterministically for that seed. Fix: a guarded SVD in the patched
  diffusionv2.py — on LinAlgError, retry with the Jacobi driver
  (gesvdj), then CPU float64. The success path is bitwise unchanged,
  so all completed arms remain valid; only runs that would have died
  take the fallback. (2) Freeze-coverage gap closed: the patched
  controller file `diffusionv2.py` — the method under test — was only
  on the GPU machine and unpinned; it now lives at the repo root and is
  covered by MANIFEST.md5.
- **2026-10-08** (mixed: one panel change, two execution fixes): (1) PANEL
  ATTRITION, decided blind to outcomes: 8 oc23-large targets
  (a0a075q0w3, a2rj53, b3eyn2, b7ie18, o76728, p21589, p33284, q9x6r4;
  474–847 aa) are excluded (`status=excluded` in targets.csv) because
  they are compute-infeasible at the frozen settings (100 samples ×
  200 steps) on 40 GB cards — they never produced a single eval CSV,
  so no diagnosis or outcome data of any kind existed when the decision
  was made; this cannot be outcome-driven attrition. Primary stratum is
  now 8 targets: a0qtt2, a6uvt1, p31133, p40131, p62495, q18a65,
  q5f9m1, q9ere7. Screened panel: 17 targets (15 oc23 + adk dev + rfah
  ms15). (2) EVAL-PATH BUG (execution-only): boltz names output
  directories after the YAML filename stem, so arms run from
  `<name>.msa.yaml` / `<name>.d32.yaml` wrote `boltz_results_<name>.msa`
  / `boltz_results_<name>.d32` while the driver looked for
  `boltz_results_<name>`. Every intervene eval and every degrade-stock
  eval failed on the path mismatch; the predictions themselves were
  valid and are recovered by re-evaluation (eval is deterministic given
  structures; the driver now skips prediction when the predictions dir
  already exists). No result was altered or lost. (3) MSA FREEZE REPAIR
  (execution-only): boltz's a3m writer appends one trailing NUL byte;
  boltz's parse_a3m crashes on it for MSAs smaller than max_seqs (the
  one larger MSA, p00558, was truncated before EOF and survived — which
  is why only its arms ran). All 17 frozen MSAs were repaired by
  stripping the single trailing NUL (content 100% intact, verified:
  exactly one NUL, at EOF, after the final newline); post-repair md5s
  logged cluster-side. freeze_msa now strips this artifact at freeze
  time and verifies, and msa_yaml refuses to feed boltz a NUL-containing
  a3m (skips the target loudly instead of falling back to the MSA
  server, which would break the no-cross-arm-MSA-drift freeze).
