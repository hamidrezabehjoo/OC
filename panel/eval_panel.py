#!/usr/bin/env python3
"""Evaluate panel predictions against two conformational references.

Panel variant of the single-target ADK/MBP evaluator (superseded; in git
history). The only difference is atom matching: predictions and references
are matched by RESIDUE NUMBER (intersection), so references with missing
density (loops, termini) or construct offsets work without any length
assertion. Everything else -- Kabsch superposition, TM-score
(d0 = 1.24*(L-15)^(1/3) - 1.8, 5 weighted iterations), basin rule
(r <= escape_rmsd and closer than the other reference) -- is identical to
the validated ADK pipeline.

References here are the renumbered files written by fetch_panel.py
(panel_data/<name>_ref1.pdb / _ref2.pdb, chain A, input-sequence numbering),
so resseq values are directly comparable with the Boltz predictions.

For fold-switch targets (foldswitch=1 in targets.csv), RMSD basins are
meaningless (the two references are different folds); pass --basin_rule tm
to assign basins from TM-scores instead:
    state1: tm1 >= 0.5 and tm1 - tm2 >= 0.15
    state2: tm2 >= 0.5 and tm2 - tm1 >= 0.15
    else:   intermediate/other
(basin column values are still written as open/closed = ref1/ref2 for
downstream compatibility with the diagnosis/analysis chain.)

Usage:
    python3 panel/eval_panel.py \
        --pred_dir results_p40131_stock_s42/boltz_results_p40131/predictions/p40131 \
        --ref1 panel_data/p40131_ref1.pdb --ref2 panel_data/p40131_ref2.pdb \
        --tag stock_s42 --out eval_p40131_stock_s42.csv
"""

import argparse
import glob
import os
import sys

import numpy as np


def read_ca_coords(path):
    """C-alpha coords of chain A keyed by residue number."""
    coords = {}
    with open(path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            if line[12:16].strip() != "CA":
                continue
            chain = line[21].strip()
            if chain and chain != "A":
                continue
            try:
                resseq = int(line[22:26])
            except ValueError:
                continue
            coords[resseq] = (float(line[30:38]), float(line[38:46]),
                              float(line[46:54]))
    if not coords:
        raise ValueError(f"no CA atoms found in {path}")
    return coords


def kabsch_superpose(mobile, ref, weights=None):
    if weights is None:
        weights = np.ones(len(ref))
    w = weights / weights.sum()
    cm = (mobile * w[:, None]).sum(axis=0)
    cr = (ref * w[:, None]).sum(axis=0)
    m = mobile - cm
    r = ref - cr
    cov = (m * w[:, None]).T @ r
    v, _, wt = np.linalg.svd(cov)
    d = np.sign(np.linalg.det(v @ wt))
    rot = v @ np.diag([1.0, 1.0, d]) @ wt
    return m @ rot + cr


def rmsd(a, b):
    return float(np.sqrt(((a - b) ** 2).sum(axis=1).mean()))


def tm_score(mobile, ref, n_iter=5):
    L = len(ref)
    d0 = max(1.24 * (L - 15) ** (1.0 / 3.0) - 1.8, 0.5)
    aln = kabsch_superpose(mobile, ref)
    for _ in range(n_iter):
        d = np.sqrt(((aln - ref) ** 2).sum(axis=1))
        w = 1.0 / (1.0 + (d / d0) ** 2)
        aln = kabsch_superpose(mobile, ref, weights=w)
    d = np.sqrt(((aln - ref) ** 2).sum(axis=1))
    return float((1.0 / (1.0 + (d / d0) ** 2)).sum() / L), rmsd(aln, ref)


def matched_arrays(pred_coords, ref_coords):
    shared = sorted(set(pred_coords) & set(ref_coords))
    if len(shared) < 20:
        return None, None, 0
    p = np.array([pred_coords[k] for k in shared])
    r = np.array([ref_coords[k] for k in shared])
    return p, r, len(shared)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pred_dir", required=True)
    ap.add_argument("--ref1", required=True, help="reference 1 (CSV column: open)")
    ap.add_argument("--ref2", required=True, help="reference 2 (CSV column: closed)")
    ap.add_argument("--tag", default="run")
    ap.add_argument("--out", default=None)
    ap.add_argument("--escape_rmsd", default="3.5",
                    help="basin RMSD threshold in A, or 'auto': min(3.5, "
                         "0.5 x ref1-ref2 RMSD) -- scaled rule for reference "
                         "pairs closer than the fixed threshold (amendment "
                         "2026-10-04)")
    ap.add_argument("--basin_rule", choices=["rmsd", "tm"], default="rmsd")
    args = ap.parse_args()

    ref1 = read_ca_coords(args.ref1)
    ref2 = read_ca_coords(args.ref2)

    if str(args.escape_rmsd).lower() == "auto":
        shared = sorted(set(ref1) & set(ref2))
        if len(shared) >= 3:
            a = np.array([ref1[k] for k in shared])
            b = np.array([ref2[k] for k in shared])
            ref_rmsd = rmsd(kabsch_superpose(a, b), b)
        else:
            ref_rmsd = float("nan")
        esc = 3.5 if not (ref_rmsd == ref_rmsd) else min(3.5, 0.5 * ref_rmsd)
        print(f"  auto escape_rmsd: ref1-ref2 RMSD {ref_rmsd:.2f} A -> threshold {esc:.2f} A")
    else:
        esc = float(args.escape_rmsd)

    if not os.path.isdir(args.pred_dir):
        sys.exit(f"prediction directory not found: {args.pred_dir}")
    files = sorted(glob.glob(
        os.path.join(args.pred_dir, "**", "*_model_*.pdb"), recursive=True))
    if not files:
        sys.exit(f"no *_model_*.pdb files under {args.pred_dir}")

    rows = []
    for path in files:
        pred = read_ca_coords(path)
        p1, r1, n1 = matched_arrays(pred, ref1)
        p2, r2, n2 = matched_arrays(pred, ref2)
        if n1 == 0 or n2 == 0:
            print(f"WARNING: {os.path.basename(path)} shares too few residues, skipped")
            continue
        tm1, ro1 = tm_score(p1, r1)
        tm2, ro2 = tm_score(p2, r2)
        if args.basin_rule == "tm":
            if tm1 >= 0.5 and tm1 - tm2 >= 0.15:
                basin = "open"
            elif tm2 >= 0.5 and tm2 - tm1 >= 0.15:
                basin = "closed"
            else:
                basin = "intermediate/other"
        else:
            if ro1 <= esc and ro1 < ro2:
                basin = "open"
            elif ro2 <= esc and ro2 < ro1:
                basin = "closed"
            else:
                basin = "intermediate/other"
        rows.append(dict(model=os.path.basename(path),
                         rmsd_open=ro1, tm_open=tm1,
                         rmsd_closed=ro2, tm_closed=tm2,
                         n1=n1, n2=n2, basin=basin))

    out = args.out or f"eval_{args.tag}.csv"
    with open(out, "w") as f:
        f.write("tag,model,rmsd_open,tm_open,rmsd_closed,tm_closed,"
                "n_matched_open,n_matched_closed,basin\n")
        for r in rows:
            f.write(f"{args.tag},{r['model']},{r['rmsd_open']:.3f},{r['tm_open']:.4f},"
                    f"{r['rmsd_closed']:.3f},{r['tm_closed']:.4f},"
                    f"{r['n1']},{r['n2']},{r['basin']}\n")

    basins = {}
    for r in rows:
        basins[r["basin"]] = basins.get(r["basin"], 0) + 1
    bo = min(rows, key=lambda r: r["rmsd_open"])
    bc = min(rows, key=lambda r: r["rmsd_closed"])
    print(f"\n=== {args.tag}: {len(rows)} models ===")
    print(f"  basin counts: {basins}")
    print(f"  best to ref1: {bo['model']}  RMSD {bo['rmsd_open']:.2f} A, TM {bo['tm_open']:.3f}")
    print(f"  best to ref2: {bc['model']}  RMSD {bc['rmsd_closed']:.2f} A, TM {bc['tm_closed']:.3f}")
    print(f"  mean RMSD to ref1 {np.mean([r['rmsd_open'] for r in rows]):.2f} A | "
          f"to ref2 {np.mean([r['rmsd_closed'] for r in rows]):.2f} A")
    print(f"  wrote {out}")


if __name__ == "__main__":
    main()
