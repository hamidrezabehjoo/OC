# Pre-written outcome paragraphs (2026-10-05(e))

Drafted before the screen, so the interpretation cannot be fitted to the
result. Fill the bracketed numbers from `panel/panel_summary.csv` /
stdout of `analyze_panel.py`. Exactly one of cases A–C applies.

## Case A — pooled diagnosed set n < 5 (inconclusive)

> The prospective screen diagnosed [n] targets in the primary stratum
> (oc23-large) as intervenable ([n_damp] damp, [n_boost] boost). Under the
> pre-registered power rule (one-sided exact signed-rank tests require
> n ≥ 5 to reach p < 0.05), the panel test is inconclusive: [describe
> shifts and per-target successes descriptively]. We make no claim of
> success or failure for the sign law on this panel; the per-target
> results, the healthy-control no-harm outcomes, and the
> MSA-depth-engineered boost arm are reported as descriptive observations
> to power future panels. [If applicable:] The pre-registered geometry
> diagnostic D1 [did / did not] flag a geometry confound
> (group p = [p], reachability: [k] targets with an unreachable state).

## Case B — efficacy passes, sign-law fails (generic broadening)

> Among the [n] diagnosed primary-stratum targets, controlled sampling
> shifted occupancy in the diagnosed direction (pooled one-sided
> signed-rank p = [p1], median shift [m]). However, the pre-registered
> sign-law test — the paired per-target comparison of each diagnosed arm
> against its contra arm — did not reach significance (p = [p2]), and
> contra arms met the success criterion in [kc]/[nct] cases. The
> occupancy gains are therefore consistent with generic broadening of
> the ensemble rather than with the predicted sign law: reducing (or
> amplifying) inference-time noise helped regardless of direction. We
> therefore do not claim a validated diagnosis rule; the controller
> remains useful as a generic ensemble-broadening knob, and the
> diagnosis–intervention link requires a design with sharper separation
> between arms.

## Case C — both pass

> Among the [n] diagnosed primary-stratum targets, controlled sampling
> shifted occupancy in the diagnosed direction (pooled one-sided
> signed-rank p = [p1], median shift [m], [k]/[n] per-target successes),
> and the sign-law test confirmed that the effect is direction-specific:
> diagnosed arms outperformed their own contra arms (paired one-sided
> signed-rank p = [p2]), with a contra false-positive rate of [kc]/[nct].
> [If D1 fired:] The geometry diagnostic flagged a correlation between
> diagnosis and basin geometry; the co-primary TM-rule analysis was
> Holm-corrected and [passed / failed] (p_TM = [p3]). The engineered
> boost arm (MSA depth 32) validated the manipulation in [k]/[n]
> healthy controls. Together these results support the claim that a
> reference-informed, two-state stock-ensemble diagnosis predicts the
> productive sign of inference-time noise control.

## Limitations (all cases)

> Our endpoint measures recovery of the deposited apo reference state:
> `expected_state` is frozen as the first reference (ref1) for every
> target, so a successful intervention means the ensemble re-occupies
> the deposited apo conformation — not necessarily the dominant state in
> solution, which for some targets is shifted by ligands, oligomerization,
> or crystallization conditions. The diagnosis is reference-informed
> under a two-state hypothesis, and targets where neither deposited
> state is reachable by the stock ensemble are flagged as geometric
> artifacts rather than diagnoses. For the two NMR-referenced targets
> (q9ere7, p62495) we use a medoid model rather than AFsample2's model
> choice, so numerical comparisons with that benchmark on those two
> targets are not like-for-like. Per-target inference rests on five
> seeds (per-target null rate 3.1%); the panel-level sign-law claim
> rests on the paired signed-rank test across targets, not on pooled
> within-target p-values.
