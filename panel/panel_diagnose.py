#!/usr/bin/env python3
"""Panel stock-ensemble diagnosis: frozen decision table, one target at a time.

Implements the same pre-registered decision table as the single-target
ADK/MBP pipeline (superseded; in git history), including the
ceiling-clause amendment (2026-10-01), plus one panel-specific extension:
for fold-switch targets (--basin_rule tm) basin fractions are re-derived
from the TM columns, because RMSD basins are meaningless when the two
references are different folds. All thresholds are identical to the
frozen protocol:

    intermediate fraction >= 0.50                          -> damp  (alpha 0.5)
    one basin share >= 0.50 AND mean pairwise Ca RMSD < 2A:
        dominant state is the expected state and >= 0.95   -> ceiling (no run)
        otherwise                                          -> boost (alpha 2.5)
    otherwise                                              -> healthy (no run)

Expected state: --expected open|closed (default open = ref1 = open for
OC23 apo periplasmic-binding proteins, alpha-hairpin for free RfaH,
inward for apo transporters).

Writes a diagnosis file (VERDICT/ALPHA lines, same format as MBP) and
appends a summary row to panel/panel_diagnosis.csv.

Usage:
    python3 panel/panel_diagnose.py --name p40131 \
        --csvs eval_p40131_stock_s42.csv \
        --pred_dirs results_p40131_stock_s42/boltz_results_p40131/predictions/p40131 \
        --basin_rule rmsd --expected open
"""

import argparse
import csv
import glob
import os
import sys

import numpy as np

# ---- frozen thresholds (identical to the frozen single-target table) ----
INTERMEDIATE_THRESHOLD = 0.50
DOMINANT_THRESHOLD = 0.50
PAIRWISE_RMSD_THRESHOLD = 2.0
CEILING_THRESHOLD = 0.95          # amendment 2026-10-01
BORDERLINE = 0.030                # amendment 2026-10-05b (audit round 2):
                                  # widened to >= 1 binomial SE at n=300, p=0.5
                                  # (SE = 2.9 pt). A point estimate 2.5 pt from
                                  # the threshold flips verdict ~19% of the
                                  # time; at 3.0 pt the guard covers 1 SE.
DAMP_ALPHA, BOOST_ALPHA = 0.5, 2.5
TM_BASIN, TM_MARGIN = 0.5, 0.15   # fold-switch basin rule (pre-registered)


def read_eval_csvs(paths, basin_rule):
    basins = []
    for p in paths:
        with open(p) as f:
            for r in csv.DictReader(f):
                if basin_rule == "tm":
                    t1, t2 = float(r["tm_open"]), float(r["tm_closed"])
                    if t1 >= TM_BASIN and t1 - t2 >= TM_MARGIN:
                        basins.append("open")
                    elif t2 >= TM_BASIN and t2 - t1 >= TM_MARGIN:
                        basins.append("closed")
                    else:
                        basins.append("intermediate/other")
                else:
                    basins.append(r["basin"])
    return basins


def read_ca(path):
    coords = {}
    with open(path) as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                ch = line[21].strip()
                if ch and ch != "A":
                    continue
                try:
                    rs = int(line[22:26])
                except ValueError:
                    continue
                coords[rs] = (float(line[30:38]), float(line[38:46]),
                              float(line[46:54]))
    idx = sorted(coords)
    return np.array([coords[i] for i in idx])


def kabsch_rmsd(a, b):
    ac, bc = a - a.mean(0), b - b.mean(0)
    v, _, wt = np.linalg.svd(ac.T @ bc)
    d = np.sign(np.linalg.det(v @ wt))
    rot = v @ np.diag([1.0, 1.0, d]) @ wt
    return float(np.sqrt(((ac @ rot - bc) ** 2).sum(1).mean()))


def mean_pairwise_rmsd(pred_dirs, max_models=100):
    xyz = []
    for d in pred_dirs:
        for p in sorted(glob.glob(os.path.join(d, "**", "*_model_*.pdb"),
                                  recursive=True))[:max_models]:
            xyz.append(read_ca(p))
    n = len(xyz)
    if n < 2:
        return float("nan"), n
    L = min(len(x) for x in xyz)
    vals = []
    for i in range(n):
        for j in range(i + 1, n):
            vals.append(kabsch_rmsd(xyz[i][:L], xyz[j][:L]))
    return float(np.mean(vals)), n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True)
    ap.add_argument("--csvs", nargs="+", required=True)
    ap.add_argument("--pred_dirs", nargs="+", required=True)
    ap.add_argument("--basin_rule", choices=["rmsd", "tm"], default="rmsd")
    ap.add_argument("--expected", choices=["open", "closed"], default="open")
    ap.add_argument("--out", default=None)
    ap.add_argument("--report_csv", default=os.path.join("panel", "panel_diagnosis.csv"))
    args = ap.parse_args()

    basins = read_eval_csvs(args.csvs, args.basin_rule)
    n = len(basins)
    fo = basins.count("open") / n
    fc = basins.count("closed") / n
    fi = basins.count("intermediate/other") / n
    fexp = fo if args.expected == "open" else fc

    # ---- diagnosis stability (amendment 2026-10-05) ----------------------
    # Two guards, both reported in the diagnosis file:
    #   (a) borderline guard: point estimate within BORDERLINE of its
    #       threshold -> indeterminate (coin-flip diagnosis at n=300)
    #   (b) per-seed agreement: the point-estimate rule is re-run on each
    #       seed's CSV alone; if the verdicts disagree -> indeterminate
    # Bootstrap CIs are reported for transparency but do not gate: at any
    # practical sample size a 95% CI straddles 0.50 for exactly the
    # near-threshold cases the panel cares about, and the error a borderline
    # guard must prevent is a coin-flip *decision*, not a wide interval.
    rng = np.random.default_rng(0)
    b = rng.choice(basins, size=(2000, n), replace=True)
    fi_ci = np.quantile((b == "intermediate/other").mean(1), [0.025, 0.975])
    dom_ci = np.quantile(np.maximum((b == "open").mean(1), (b == "closed").mean(1)),
                         [0.025, 0.975])

    def point_verdict(bb):
        _fi = bb.count("intermediate/other") / len(bb)
        _dom = max(bb.count("open"), bb.count("closed")) / len(bb)
        if _fi >= INTERMEDIATE_THRESHOLD:
            return "damp"
        return "collapsed?" if _dom >= DOMINANT_THRESHOLD else "healthy"

    per_seed = [point_verdict(read_eval_csvs([c], args.basin_rule)) for c in args.csvs]
    unanimous = len(set(per_seed)) == 1

    lines = [f"panel diagnosis for {args.name}",
             f"stock ensemble: {n} models ({len(args.csvs)} seed set(s))",
             f"  ref1/open fraction        {fo:.3f}",
             f"  ref2/closed fraction      {fc:.3f}",
             f"  intermediate fraction     {fi:.3f}  (95% CI {fi_ci[0]:.3f}-{fi_ci[1]:.3f})",
             f"  dominant-basin share      {max(fo,fc):.3f}  (95% CI {dom_ci[0]:.3f}-{dom_ci[1]:.3f})",
             f"  per-seed verdicts: {per_seed}",
             f"  basin rule {args.basin_rule}, expected state {args.expected}", ""]

    borderline = abs(fi - INTERMEDIATE_THRESHOLD) <= BORDERLINE or \
                 (fi < INTERMEDIATE_THRESHOLD and
                  abs(max(fo, fc) - DOMINANT_THRESHOLD) <= BORDERLINE)
    if borderline or not unanimous:
        why = []
        if borderline:
            why.append(f"point estimate within {BORDERLINE} of threshold")
        if not unanimous:
            why.append(f"per-seed verdicts disagree {per_seed}")
        lines.append("rule: " + " AND ".join(why) + " -> indeterminate")
        _finish(args, lines, n, fo, fc, fi, "indeterminate", "")
        return

    verdict, alpha = "healthy", ""
    if fi >= INTERMEDIATE_THRESHOLD:
        verdict, alpha = "damp", DAMP_ALPHA
        lines.append(f"rule: intermediate {fi:.3f} >= {INTERMEDIATE_THRESHOLD} -> over-diffuse")
    else:
        dom = max(fo, fc)
        if dom >= DOMINANT_THRESHOLD:
            mp, nm = mean_pairwise_rmsd(args.pred_dirs)
            lines.append(f"rule: dominant share {dom:.3f} >= {DOMINANT_THRESHOLD}; "
                         f"mean pairwise Ca RMSD {mp:.2f} A over {nm} models")
            if mp < PAIRWISE_RMSD_THRESHOLD:
                if fexp >= CEILING_THRESHOLD:
                    verdict = "ceiling"
                    lines.append(f"rule: expected-state fraction {fexp:.3f} >= "
                                 f"{CEILING_THRESHOLD} -> ceiling, no actionable mismatch")
                else:
                    verdict, alpha = "boost", BOOST_ALPHA
                    lines.append("rule: collapsed ensemble, expected state not at ceiling "
                                 "-> boost")
            else:
                lines.append("rule: dominant basin but ensemble not collapsed -> healthy")
        else:
            lines.append("rule: no branch triggered -> healthy")

    lines += ["", f"VERDICT {verdict}" + (f" ALPHA {alpha}" if alpha else "")]
    _finish(args, lines, n, fo, fc, fi, verdict, alpha)


def _finish(args, lines, n, fo, fc, fi, verdict, alpha):
    out = args.out or f"diagnosis_{args.name}.txt"
    with open(out, "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))

    # panel_diagnosis.csv: one row per target, latest diagnosis
    fields = ["name", "n", "frac_ref1", "frac_ref2", "frac_intermediate",
              "basin_rule", "expected", "verdict", "alpha"]
    old = {}
    if os.path.exists(args.report_csv):
        for r in csv.DictReader(open(args.report_csv)):
            old[r["name"]] = r
    old[args.name] = dict(zip(fields, [args.name, n, f"{fo:.3f}", f"{fc:.3f}",
                                       f"{fi:.3f}", args.basin_rule, args.expected,
                                       verdict, alpha]))
    with open(args.report_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for _, r in old.items():
            w.writerow(r)
    print(f"updated {args.report_csv}")


if __name__ == "__main__":
    main()
