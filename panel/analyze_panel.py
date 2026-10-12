#!/usr/bin/env python3
"""Panel-level analysis of the prospective screen (protocol v2.1).

Per-target evaluation (frozen success criteria, protocol section 5):

  damp  ctrl arm : expected-state occupancy increases vs stock,
                   Fisher two-sided p < 0.05, direction in 3/3 seeds.
  boost ctrl arm : minority-basin occupancy increases vs stock,
                   Fisher p < 0.05, 3/3 seeds.
  contra arm     : sign-law control. The contra-diagnosed arm (boost where
                   the rule says damp, and vice versa) is evaluated on the
                   same metric; the sign law predicts NO improvement.
  healthy control: no-harm criterion: expected-state occupancy does not
                   decrease under either arm (one-sided Fisher p >= 0.05
                   for a decrease; i.e. we fail to detect harm).

Strata (pre-registered, amendment 2026-10-05):
  dev       - ADK and MBP: development cases, reported separately, never
              in the primary tally (the diagnosis rule and endpoint were
              validated/calibrated on them)
  oc23-large- OC23 targets with reference-pair Ca RMSD >= 4.0 A  (PRIMARY)
  oc23-small- OC23 targets below 4.0 A (scaled basin threshold; secondary)
  ms15      - exploratory stratum (hand-picked rows)
  tp16      - optional membrane-transporter stratum

Panel-level endpoints:
  primary   - success fraction among intervened PRIMARY-stratum targets,
              Wilson 95% CI; sign-law contrast: fraction of contra arms
              showing the same-direction improvement
  secondary - per-stratum success, diagnosis-bin distribution, median
              occupancy shift, healthy-control no-harm rate
  power rule: a diagnosis bin with < 4 members is reported as
              inconclusive, not failed.

Usage:
    python3 panel/analyze_panel.py --targets panel/targets.csv
"""

import argparse
import csv
import glob
import math
import os

TM_BASIN, TM_MARGIN = 0.5, 0.15
LARGE_MOTION_CUTOFF = 4.0   # reference-pair Ca RMSD, amendment 2026-10-05
DEV_CASES = {"adk", "mbp"}
MIN_BIN = 5   # amendment 2026-10-05(c): one-sided exact Wilcoxon cannot
              # reach p < 0.05 below n = 5 (min p = 0.0625 at n = 4)


def logchoose(n, k):
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def fisher_two_sided(a, b, c, d):
    n = a + b + c + d
    r1, c1 = a + b, a + c
    lo, hi = max(0, r1 + c1 - n), min(r1, c1)

    def prob(x):
        return math.exp(logchoose(r1, x) + logchoose(n - r1, c1 - x)
                        - logchoose(n, c1))

    p_obs = prob(a)
    return min(1.0, sum(prob(x) for x in range(lo, hi + 1)
                        if prob(x) <= p_obs * (1 + 1e-9)))


def fisher_greater(a, b, c, d):
    """One-sided P(X >= a) for [[a,b],[c,d]] (ctrl enriched in row 1)."""
    n = a + b + c + d
    r1, c1 = a + b, a + c
    lo, hi = max(0, r1 + c1 - n), min(r1, c1)
    return min(1.0, sum(math.exp(logchoose(r1, x) + logchoose(n - r1, c1 - x)
                                 - logchoose(n, c1))
                        for x in range(a, hi + 1)))


def wilson(k, n, z=1.959964):
    if n == 0:
        return float("nan"), float("nan"), float("nan")
    p = k / n
    den = 1 + z * z / n
    cen = (p + z * z / (2 * n)) / den
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return p, cen - half, cen + half


def basins_from_csv(path, foldswitch):
    out = []
    with open(path) as f:
        for r in csv.DictReader(f):
            if foldswitch:
                t1, t2 = float(r["tm_open"]), float(r["tm_closed"])
                if t1 >= TM_BASIN and t1 - t2 >= TM_MARGIN:
                    out.append("open")
                elif t2 >= TM_BASIN and t2 - t1 >= TM_MARGIN:
                    out.append("closed")
                else:
                    out.append("intermediate/other")
            else:
                out.append(r["basin"])
    return out


def frac(b, key):
    return b.count(key) / len(b) if b else float("nan")


def min_rmsd_to_refs(stock_csvs):
    """Best-of-stock RMSD to each reference (reachability, amendment
    2026-10-05(c)): a state the stock ensemble never gets within the basin
    radius of cannot host a diagnosable basin."""
    m1 = m2 = float("inf")
    for p in stock_csvs:
        with open(p) as f:
            for r in csv.DictReader(f):
                try:
                    m1 = min(m1, float(r["rmsd_open"]))
                    m2 = min(m2, float(r["rmsd_closed"]))
                except (ValueError, KeyError):
                    pass
    return m1, m2


def load_arm(name, arm, seeds, foldswitch):
    per_seed = []
    for s in seeds:
        p = f"eval_{name}_{arm}_s{s}.csv"
        if os.path.exists(p):
            per_seed.append((s, basins_from_csv(p, foldswitch)))
    return per_seed


def seed_bootstrap_ci(stock, per_seed, key, n_boot=2000, rng_seed=0):
    """Seed-level (block) bootstrap CI of the occupancy shift.
    Resamples SEEDS with replacement; each replicate pools 100 models per
    resampled seed. This is the honest uncertainty under shared trunk/MSA
    (amendment 2026-10-05(c), pseudo-replication audit): models within a
    seed are not independent, so the pooled Fisher p is descriptive only."""
    import random
    rng = random.Random(rng_seed)
    m_s = frac(stock, key)
    seeds = [bb for _, bb in per_seed]
    shifts = []
    for _ in range(n_boot):
        draw = [seeds[rng.randrange(len(seeds))] for _ in seeds]
        pooled = [b for bb in draw for b in bb]
        shifts.append(frac(pooled, key) - m_s)
    shifts.sort()
    lo = shifts[int(0.025 * n_boot)]
    hi = shifts[int(0.975 * n_boot)]
    return lo, hi


def eval_arm(stock, per_seed, key, direction):
    """direction: 'up' (occupancy of key increases) or 'down'.
    Success (amendment 2026-10-05(c)): direction in 3/3 seeds AND the
    seed-level bootstrap 95% CI excludes 0. The pooled Fisher p is
    reported as descriptive only (within-seed predictions share trunk/MSA
    and are not independent)."""
    if not per_seed:
        return None
    pooled = [b for _, bb in per_seed for b in bb]
    m_s, m_c = frac(stock, key), frac(pooled, key)
    a, bb_ = pooled.count(key), len(pooled) - pooled.count(key)
    c, d = stock.count(key), len(stock) - stock.count(key)
    p_two = fisher_two_sided(a, bb_, c, d)
    ok_dir = (m_c > m_s) if direction == "up" else (m_c < m_s)
    seed_shifts = [(frac(x, key) - m_s) if direction == "up" else (m_s - frac(x, key))
                   for _, x in per_seed]
    n_ok = sum(s > 0 for s in seed_shifts)
    blo, bhi = seed_bootstrap_ci(stock, per_seed, key)
    if direction == "down":
        blo, bhi = -bhi, -blo
    ci_ok = blo > 0 if direction == "up" else bhi < 0
    return dict(metric_stock=m_s, metric_ctrl=m_c, shift=m_c - m_s,
                p=p_two, seeds_ok=f"{n_ok}/{len(per_seed)}",
                seed_shifts=seed_shifts, boot_lo=blo, boot_hi=bhi,
                success=ok_dir and ci_ok and n_ok == len(per_seed)
                and len(per_seed) >= 3)


def spearman(x, y):
    """Spearman rank correlation (average ranks for ties)."""
    n = len(x)
    if n < 3:
        return float("nan"), n

    def rank(v):
        order = sorted(range(n), key=lambda i: v[i])
        r = [0.0] * n
        i = 0
        while i < n:
            j = i
            while j + 1 < n and v[order[j + 1]] == v[order[i]]:
                j += 1
            for k in range(i, j + 1):
                r[order[k]] = (i + j) / 2.0 + 1.0
            i = j + 1
        return r

    rx, ry = rank(x), rank(y)
    mx, my = sum(rx) / n, sum(ry) / n
    cov = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    vx = sum((a - mx) ** 2 for a in rx)
    vy = sum((b - my) ** 2 for b in ry)
    return (cov / math.sqrt(vx * vy)) if vx > 0 and vy > 0 else float("nan"), n


def wilcoxon_signed_rank_exact(diffs, alternative="two-sided"):
    """Exact Wilcoxon signed-rank p by full enumeration (n <= 20).
    Zeros dropped; ties get average ranks. alternative: 'two-sided' or
    'greater' (H1: median > 0). Pre-registered one-sided for the panel
    primary test (amendment 2026-10-05c): the sign law is directional, and
    two-sided needs n >= 6 to reach 0.05 (min p: 0.125 @ n=4, 0.0625 @
    n=5); one-sided reaches 0.03125 at n = 5."""
    d = [x for x in diffs if x != 0]
    n = len(d)
    if n < 3:
        return float("nan"), n
    order = sorted(range(n), key=lambda i: abs(d[i]))
    ranks = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j + 1 < n and abs(abs(d[order[j + 1]]) - abs(d[order[i]])) < 1e-12:
            j += 1
        for k in range(i, j + 1):
            ranks[order[k]] = (i + j) / 2.0 + 1.0
        i = j + 1
    w_obs = sum(r for r, x in zip(ranks, d) if x > 0)
    tot = sum(ranks)
    # enumerate all 2^n sign assignments over ranks
    cnt, extreme = 0, 0
    for mask in range(1 << n):
        w = sum(ranks[i] for i in range(n) if mask >> i & 1)
        cnt += 1
        if alternative == "greater":
            if w >= w_obs - 1e-9:
                extreme += 1
        elif abs(w - tot / 2.0) >= abs(w_obs - tot / 2.0) - 1e-9:
            extreme += 1
    return extreme / cnt, n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", required=True)
    ap.add_argument("--seeds", nargs="+", default=["42", "1", "2"],
                    help="stock screen seeds (diagnosis)")
    ap.add_argument("--ctrl_seeds", nargs="+",
                    default=["42", "1", "2", "7", "13"],
                    help="controlled/contra arm seeds (amendment "
                         "2026-10-05(e): 5 seeds; at 3 the bootstrap is "
                         "degenerate and the null rate is 12.5%%)")
    ap.add_argument("--diag", default=os.path.join("panel", "panel_diagnosis.csv"))
    ap.add_argument("--fetch_report", default=os.path.join("panel", "fetch_report.csv"))
    ap.add_argument("--out", default=os.path.join("panel", "panel_summary.csv"))
    args = ap.parse_args()

    meta = {r["name"]: r for r in csv.DictReader(open(args.targets))}
    diag = {r["name"]: r for r in csv.DictReader(open(args.diag))}
    refrmsd = {}
    if os.path.exists(args.fetch_report):
        for r in csv.DictReader(open(args.fetch_report)):
            try:
                refrmsd[r["name"]] = float(r["ref_rmsd"])
            except (ValueError, KeyError):
                pass

    def stratum(name):
        if name in DEV_CASES:
            return "dev"
        st = meta.get(name, {}).get("set", "?")
        if st == "oc23":
            rr = refrmsd.get(name, float("nan"))
            if rr == rr and rr < LARGE_MOTION_CUTOFF:
                return "oc23-small"
            return "oc23-large"
        return st

    rows_out = []
    for name, dg in diag.items():
        foldswitch = meta.get(name, {}).get("foldswitch", "0") == "1"
        verdict = dg["verdict"]
        stock_csvs = sorted(glob.glob(f"eval_{name}_stock_s*.csv"))
        if not stock_csvs:
            continue
        stock = []
        for p in stock_csvs:
            stock += basins_from_csv(p, foldswitch)
        s1, s2 = frac(stock, "open"), frac(stock, "closed")
        exp_key = "open" if dg["expected"] == "open" else "closed"
        minority = "open" if s1 <= s2 else "closed"

        rec = dict(name=name, stratum=stratum(name), verdict=verdict,
                   alpha=dg["alpha"], n_stock=len(stock),
                   stock_ref1=f"{s1:.3f}", stock_ref2=f"{s2:.3f}",
                   stock_inter=f"{frac(stock,'intermediate/other'):.3f}")

        # reachability annotation (amendment 2026-10-05(c))
        rr = refrmsd.get(name, float("nan"))
        if rr == rr:
            radius = min(3.5, 0.5 * rr)
            mo, mc = min_rmsd_to_refs(stock_csvs)
            flags = []
            if mo > radius:
                flags.append(f"ref1 unreachable (best {mo:.2f} > {radius:.2f})")
            if mc > radius:
                flags.append(f"ref2 unreachable (best {mc:.2f} > {radius:.2f})")
            rec["reach"] = "; ".join(flags) if flags else "both reachable"
            if len(flags) == 2 and verdict in ("damp", "boost"):
                rec["reach"] += " <-- BOTH unreachable: diagnosis is geometric, treat as indeterminate"

        if verdict in ("damp", "boost"):
            key = exp_key if verdict == "damp" else minority
            ctrl = eval_arm(stock, load_arm(name, "ctrl", args.ctrl_seeds, foldswitch),
                            key, "up")
            contra = eval_arm(stock, load_arm(name, "contra", args.ctrl_seeds, foldswitch),
                              key, "up")
            # TM-rule sensitivity shifts (co-primary under Holm only if
            # diagnostic D1 fires; amendment 2026-10-05(e)). Standard
            # length-dependent d0 — the fixed d0 = 3.5 A convention was
            # used ONLY to match AFsample2's table (their TMalign -d 3.5),
            # never for the sensitivity re-diagnosis.
            stock_tm = []
            for p in stock_csvs:
                stock_tm += basins_from_csv(p, True)
            s1t, s2t = frac(stock_tm, "open"), frac(stock_tm, "closed")
            minority_tm = "open" if s1t <= s2t else "closed"
            key_tm = exp_key if verdict == "damp" else minority_tm
            ctrl_tm = eval_arm(stock_tm, load_arm(name, "ctrl", args.ctrl_seeds, True),
                               key_tm, "up")
            contra_tm = eval_arm(stock_tm, load_arm(name, "contra", args.ctrl_seeds, True),
                                 key_tm, "up")
            if ctrl is None:
                rec.update(success="", detail="controlled arm not run yet")
            else:
                rec.update(n_ctrl="", metric_stock=f"{ctrl['metric_stock']:.3f}",
                           metric_ctrl=f"{ctrl['metric_ctrl']:.3f}",
                           shift=f"{ctrl['shift']:+.3f}",
                           fisher_p=f"{ctrl['p']:.2e}", seeds_ok=ctrl["seeds_ok"],
                           shift_seeds=";".join(f"{s:+.3f}" for s in ctrl["seed_shifts"]),
                           boot_ci=f"[{ctrl['boot_lo']:+.3f},{ctrl['boot_hi']:+.3f}]",
                           success="yes" if ctrl["success"] else "no",
                           contra_shift=(f"{contra['shift']:+.3f}" if contra else ""),
                           contra_p=(f"{contra['p']:.2e}" if contra else ""),
                           contra_success=("" if contra is None
                                           else "yes" if contra["success"] else "no"),
                           detail="")
                rec["_ctrl_shift"] = ctrl["shift"]
                if contra is not None:
                    rec["_contra_shift"] = contra["shift"]
                if ctrl_tm is not None:
                    rec["_ctrl_shift_tm"] = ctrl_tm["shift"]
                if contra_tm is not None:
                    rec["_contra_shift_tm"] = contra_tm["shift"]
        elif verdict == "healthy":
            # no-harm control (if arms were run): expected state must not drop
            harm = []
            for arm in ("ctrl", "contra"):
                per = load_arm(name, arm, args.ctrl_seeds, foldswitch)
                if not per:
                    continue
                pooled = [b for _, bb in per for b in bb]
                a, bb_ = pooled.count(exp_key), len(pooled) - pooled.count(exp_key)
                c, d = stock.count(exp_key), len(stock) - stock.count(exp_key)
                # one-sided: stock enriched vs arm => arm harmed expected state
                p_harm = fisher_greater(c, d, a, bb_)
                harm.append(f"{arm}: dExp={frac(pooled,exp_key)-frac(stock,exp_key):+.3f} "
                            f"p_harm={p_harm:.2e}")
            rec.update(success="", detail="; ".join(harm) if harm
                       else "no intervention (protocol)")
        else:
            rec.update(success="", detail="no intervention (protocol)")
        rows_out.append(rec)

    fields = ["name", "stratum", "verdict", "alpha", "n_stock", "stock_ref1",
              "stock_ref2", "stock_inter", "n_ctrl", "metric_stock",
              "metric_ctrl", "shift", "shift_seeds", "boot_ci", "fisher_p",
              "seeds_ok", "success", "contra_shift", "contra_p",
              "contra_success", "reach", "detail"]
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for r in rows_out:
            w.writerow(r)

    # ---------------- panel-level report ----------------
    bins = {}
    for r in rows_out:
        bins[r["verdict"]] = bins.get(r["verdict"], 0) + 1
    done = [r for r in rows_out if r.get("success") in ("yes", "no")]

    # ---- pre-registered diagnostic D1 (amendment 2026-10-05c): does the
    # diagnosis track geometry rather than physics? The basin radius is
    # capped (min(3.5, 0.5 x ref-sep)), so radius/sep shrinks for wide
    # hinges (ADK 0.49 vs q9ere7 0.17) and borderline models fall outside
    # both basins -> mechanical 'intermediate'. Check Spearman correlation
    # of the stock intermediate fraction with radius/sep and with length.
    d1 = []
    for r in rows_out:
        rr = refrmsd.get(r["name"])
        try:
            L = float(meta.get(r["name"], {}).get("length", "nan"))
        except ValueError:
            L = float("nan")
        if rr and rr == rr and L == L and r.get("stock_inter"):
            radius = min(3.5, 0.5 * rr)
            d1.append((r["name"], radius / rr, L, float(r["stock_inter"]), rr))
    if len(d1) >= 6:
        rho_rs, n1 = spearman([d[1] for d in d1], [d[3] for d in d1])
        rho_L, n2 = spearman([d[2] for d in d1], [d[3] for d in d1])
        flag = " <-- GEOMETRY CONFOUND (see protocol section 6)" if abs(rho_rs) >= 0.5 else ""
        print(f"\nD1 geometry check: Spearman(intermediate fraction, "
              f"radius/sep) = {rho_rs:+.2f}, vs length = {rho_L:+.2f} "
              f"(n={n1}){flag}")
        # group comparison (amendment 2026-10-05(c)): radius/sep is pinned
        # at 0.5 for every target with separation <= 7 A (~19 of 25), so the
        # Spearman has little power; compare the two separation groups
        # directly with an exact two-sided permutation test on medians.
        lo = sorted(d[3] for d in d1 if d[4] <= 7.0)
        hi = sorted(d[3] for d in d1 if d[4] > 7.0)
        if lo and hi:
            def med(v): return v[len(v) // 2] if len(v) % 2 else 0.5 * (v[len(v) // 2 - 1] + v[len(v) // 2])
            obs = med(hi) - med(lo)
            allv = lo + hi
            k = len(hi)
            from itertools import combinations
            n_comb = 0
            extr = 0
            idx = range(len(allv))
            if len(allv) <= 25:
                for comb in combinations(idx, k):
                    s = sum(allv[i] for i in comb)
                    mh = s / k
                    ml = (sum(allv) - s) / len(lo)
                    n_comb += 1
                    if abs(mh - ml) >= abs(obs) - 1e-12:
                        extr += 1
                p_perm = extr / n_comb
                print(f"D1 group check: sep>7 A (n={len(hi)}) median intermediate "
                      f"{med(hi):.3f} vs sep<=7 A (n={len(lo)}) {med(lo):.3f}, "
                      f"exact permutation p={p_perm:.4f}"
                      + (" <-- GEOMETRY CONFOUND" if p_perm < 0.05 and abs(obs) > 0.1 else ""))
        # reachability tally
        unreach = [r["name"] for r in rows_out
                   if r.get("reach", "").startswith("ref") or "BOTH" in r.get("reach", "")]
        if unreach:
            print(f"D1 reachability: {len(unreach)} targets with an "
                  f"unreachable reference: {', '.join(sorted(unreach))}")

    print("\n================ PANEL SUMMARY ================")
    print(f"screened targets: {len(rows_out)}")
    print(f"diagnosis bins: {bins}")

    # PANEL TESTS (amendment 2026-10-05c, roles fixed):
    #   EFFICACY CLAIM  = pooled one-sided Wilcoxon of ctrl-arm shifts vs 0,
    #                     signed in the diagnosed direction. (Can pass from
    #                     generic broadening alone — not the sign-law claim.)
    #   SIGN-LAW CLAIM  = paired ctrl-vs-contra one-sided Wilcoxon: the
    #                     diagnosed arm must beat its contra arm per target.
    # damp/boost branch fractions are descriptive only (bins < 5 cannot
    # reach significance even one-sided; see protocol section 6).
    prim_diag = [r for r in done if r["stratum"] == "oc23-large"
                 and r["verdict"] in ("damp", "boost") and "_ctrl_shift" in r]
    if len(prim_diag) >= 3:
        shifts_d = [r["_ctrl_shift"] for r in prim_diag]
        p1, n1 = wilcoxon_signed_rank_exact(shifts_d, alternative="greater")
        k = sum(r["success"] == "yes" for r in prim_diag)
        wp, wl, wh = wilson(k, len(prim_diag))
        med = sorted(shifts_d)[len(shifts_d) // 2]
        print(f"EFFICACY (oc23-large, pooled diagnosed): {k}/{len(prim_diag)} "
              f"per-target successes ({wp:.2f}, Wilson {wl:.2f}-{wh:.2f}); "
              f"median signed shift {med:+.3f}; "
              f"one-sided exact Wilcoxon p={p1:.4f} (n={n1}) "
              f"-> {'SIGNIFICANT' if p1 == p1 and p1 < 0.05 else 'n.s.'}")
        paired = [(r["_ctrl_shift"], r["_contra_shift"]) for r in prim_diag
                  if "_contra_shift" in r]
        if len(paired) >= 3:
            pc1, nc1 = wilcoxon_signed_rank_exact([c - k_ for c, k_ in paired],
                                                  alternative="greater")
            kc = sum(r.get("contra_success") == "yes" for r in prim_diag
                     if r.get("contra_success") in ("yes", "no"))
            nct = sum(r.get("contra_success") in ("yes", "no") for r in prim_diag)
            print(f"SIGN-LAW (primary claim): paired ctrl-vs-contra one-sided "
                  f"Wilcoxon p={pc1:.4f} (n={nc1}) "
                  f"-> {'SIGNIFICANT' if pc1 == pc1 and pc1 < 0.05 else 'n.s.'}; "
                  f"contra false-positive rate {kc}/{nct}")
            # Holm correction across the two basin-rule analyses, needed
            # only if D1 fires (TM rule becomes co-primary; amendment
            # 2026-10-05(e)). Computed always, reported always.
            paired_tm = [(r["_ctrl_shift_tm"], r["_contra_shift_tm"])
                         for r in prim_diag
                         if "_ctrl_shift_tm" in r and "_contra_shift_tm" in r]
            if len(paired_tm) >= 3:
                p_tm, _ = wilcoxon_signed_rank_exact(
                    [c - k_ for c, k_ in paired_tm], alternative="greater")
                ps = sorted([pc1, p_tm])
                holm = ps[0] < 0.025 or (ps[0] < 0.05 and ps[1] < 0.05)
                print(f"  TM-rule paired p={p_tm:.4f}; Holm over "
                      f"{{RMSD-rule, TM-rule}} (applies if D1 flags): "
                      f"{'PASS' if holm else 'FAIL'}")
        for br in ("damp", "boost"):
            sub = [r for r in prim_diag if r["verdict"] == br]
            if sub:
                ks = sum(r["success"] == "yes" for r in sub)
                ms = sorted(r["_ctrl_shift"] for r in sub)[len(sub) // 2]
                print(f"  {br} branch (descriptive): {ks}/{len(sub)} "
                      f"successes, median shift {ms:+.3f}")
        # robustness subset (amendment 2026-10-05c): primary targets with
        # ref-pair RMSD >= 5.0 A — free of boundary-fragility worries
        # (4.0 A sits only 0.23 A above the MBP calibration point)
        xl = [r for r in prim_diag if refrmsd.get(r["name"], 0.0) >= 5.0]
        if len(xl) >= 3:
            px, nx = wilcoxon_signed_rank_exact([r["_ctrl_shift"] for r in xl],
                                                alternative="greater")
            kx = sum(r["success"] == "yes" for r in xl)
            print(f"  robustness subset (>=5 A, n={len(xl)}): {kx}/{len(xl)} "
                  f"successes, one-sided Wilcoxon p={px:.4f}")
    elif prim_diag:
        print(f"PRIMARY: only {len(prim_diag)} diagnosed primary targets "
              f"-- inconclusive (<3)")
    for st in ("oc23-large", "oc23-small", "ms15", "tp16", "dev"):
        sub = [r for r in done if r["stratum"] == st]
        if not sub:
            continue
        k = sum(r["success"] == "yes" for r in sub)
        tag = " (PRIMARY)" if st == "oc23-large" else ""
        if len(sub) < MIN_BIN and st != "dev":
            print(f"  {st}: {k}/{len(sub)} -- inconclusive (<{MIN_BIN} members){tag}")
            continue
        p, lo, hi = wilson(k, len(sub))
        print(f"  {st}: {k}/{len(sub)} = {p:.2f} (Wilson 95% CI {lo:.2f}-{hi:.2f}){tag}")
        contra = [r for r in sub if r.get("contra_success") in ("yes", "no")]
        if contra:
            kc = sum(r["contra_success"] == "yes" for r in contra)
            pc, loc, hic = wilson(kc, len(contra))
            print(f"    sign-law control: contra false-positive rate "
                  f"{kc}/{len(contra)} = {pc:.2f} (Wilson 95% CI {loc:.2f}-{hic:.2f}); "
                  f"sign law predicts 0")
        # primary comparison: paired per-target diagnosed-arm shift vs
        # contra-arm shift (audit round 2, point 4)
        pairs = [(r["_ctrl_shift"], r["_contra_shift"]) for r in sub
                 if "_ctrl_shift" in r and "_contra_shift" in r]
        if len(pairs) >= 3:
            pw, npair = wilcoxon_signed_rank_exact([c - k_ for c, k_ in pairs],
                                                   alternative="greater")
            mc = sorted(c for c, _ in pairs)[len(pairs) // 2]
            mk = sorted(k_ for _, k_ in pairs)[len(pairs) // 2]
            print(f"    paired ctrl-vs-contra: median shift {mc:+.3f} vs {mk:+.3f}, "
                  f"one-sided Wilcoxon (exact) p={pw:.4f} (n={npair})")
        elif pairs:
            print(f"    paired ctrl-vs-contra: only {len(pairs)} paired targets "
                  f"-- inconclusive (<3 pairs)")
    dev = [r for r in done if r["stratum"] == "dev"]
    if dev:
        for r in dev:
            print(f"  dev case {r['name']}: {r['verdict']} -> success={r['success']} "
                  f"(reported separately, not in primary tally)")
    shifts = sorted(float(r["shift"]) for r in done if r["success"] == "yes")
    if shifts:
        print(f"median occupancy shift (successful): {shifts[len(shifts)//2]:+.3f}")

    # ---- engineered boost: MSA-depth degradation (secondary arm,
    # amendments 2026-10-05b/c) ----
    # Manipulation check: an arm counts as valid only if the DEGRADED
    # stock (depth-32 MSA, no control) re-diagnoses as 'boost'
    # (collapsed into one basin); the alpha=2.5 arm is then evaluated
    # against the degraded stock, not the full-MSA stock. Pass rate
    # reported.
    ddiag = {}
    dpath = os.path.join("panel", "panel_diagnosis_degrade.csv")
    if os.path.exists(dpath):
        for r in csv.DictReader(open(dpath)):
            ddiag[r["name"]] = r
    deg_names = set()
    for p in glob.glob("eval_*_degrade_s*.csv"):
        deg_names.add(p[len("eval_"):p.index("_degrade_s")])
    n_valid = 0
    for name in sorted(deg_names):
        fs = meta.get(name, {}).get("foldswitch", "0") == "1"
        dstock = []
        for p in sorted(glob.glob(f"eval_{name}_degrade_stock_s*.csv")):
            dstock += basins_from_csv(p, fs)
        dv = ddiag.get(name, {}).get("verdict", "?")
        valid = dv == "boost"
        n_valid += valid
        line = (f"  engineered boost (MSA depth 32) on {name}: "
                f"degraded-stock diagnosis = {dv} "
                f"({'manipulation VALID' if valid else 'manipulation FAILED'})")
        if valid and dstock:
            pooled = []
            for p in sorted(glob.glob(f"eval_{name}_degrade_s[0-9]*.csv")):
                pooled += basins_from_csv(p, fs)
            if pooled:
                s1, s2 = frac(dstock, "open"), frac(dstock, "closed")
                minority = "open" if s1 <= s2 else "closed"
                m_s, m_d = frac(dstock, minority), frac(pooled, minority)
                a, bb_ = pooled.count(minority), len(pooled) - pooled.count(minority)
                c, d = dstock.count(minority), len(dstock) - dstock.count(minority)
                p_two = fisher_two_sided(a, bb_, c, d)
                line += (f"; minority basin {m_s:.3f} -> {m_d:.3f} "
                         f"(p={p_two:.2e}) "
                         f"{'BOOST-SUCCESS' if m_d > m_s and p_two < 0.05 else ''}")
        print(line)
    if deg_names:
        print(f"  manipulation-check pass rate: {n_valid}/{len(deg_names)}")
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
