#!/usr/bin/env python3
"""Fetch and prepare panel targets: sequences, references, boltz yamls.

For every row of targets.csv with status=ready (and not optional, unless
--include-optional), this script:

  1. obtains the input sequence (UniProt FASTA if the row has a uniprot id,
     otherwise the SEQRES of the ref1 chain from the PDB);
  2. downloads both reference structures from RCSB and extracts the
     specified chain;
  3. aligns the reference chain sequence to the input sequence
     (Needleman-Wunsch) and rewrites the reference as
     panel_data/<name>_ref1.pdb / <name>_ref2.pdb containing only aligned
     residues, chain id A, renumbered to input-sequence positions (1..L);
  4. writes yamls/<name>.yaml (same schema as adk.yaml / mbp.yaml);
  5. runs sanity checks and appends one row to fetch_report.csv.

Sanity checks (all reported, none fatal except download failure):
  - n_matched: residues of each reference aligned to the input sequence
  - shared: positions covered by BOTH references (what eval_panel.py uses)
  - ref-rmsd: Ca RMSD between the two references on shared positions
      -> WARNING if < 3.5 A (the basin threshold is degenerate there)
  - coverage warning if shared < 60% of input length

Requires only numpy + network access. No biopython.

Usage:
    python3 panel/fetch_panel.py --targets panel/targets.csv
    python3 panel/fetch_panel.py --targets panel/targets.csv --only adk,rfah
    python3 panel/fetch_panel.py --targets panel/targets.csv --include-optional
"""

import argparse
import csv
import os
import sys
import urllib.request

import numpy as np

AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M",  # selenomethionine -> methionine
}


def fetch(url, dest=None, binary=False, retries=4):
    import time
    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "panel-fetch/1.0"})
            data = urllib.request.urlopen(req, timeout=180).read()
            if not data:
                raise IOError("empty response")
            if dest:
                with open(dest, "wb") as f:
                    f.write(data)
            return data if binary else data.decode()
        except Exception as e:
            last = e
            if dest and os.path.exists(dest):
                os.remove(dest)   # never leave a partial cache file
            time.sleep(2 * (attempt + 1))
    raise last


def uniprot_sequence(uid):
    fasta = fetch(f"https://rest.uniprot.org/uniprotkb/{uid}.fasta")
    return "".join(l.strip() for l in fasta.splitlines() if not l.startswith(">"))


# Filled by parse_pdb_chain: path -> "medoid model k of n" when an NMR
# ensemble was reduced, for the fetch report.
MODEL_PICK = {}


def pick_medoid(models):
    """Pre-registered NMR reference choice (amendment 2026-10-05b): the model
    minimizing mean pairwise Ca RMSD to all other models over shared residues;
    ties -> lowest model number."""
    common = set(models[0][0])
    for res, _ in models[1:]:
        common &= set(res)
    if len(common) < 10:
        return 0
    idx = sorted(common)
    coords = [np.array([m[0][r][1] for r in idx]) for m in models]
    best_i, best_s = 0, np.inf
    for i in range(len(models)):
        s = sum(kabsch_rmsd(coords[i], coords[j])
                for j in range(len(models)) if j != i)
        s /= max(len(models) - 1, 1)
        if s < best_s - 1e-9:
            best_i, best_s = i, s
    return best_i


def parse_pdb_chain(path, chain_want):
    """Return (seqres_seq, {resseq: (resname, ca_xyz)}, atom_lines).

    NMR multi-model entries: every model is parsed and the medoid model
    (see pick_medoid) is used instead of MODEL 1."""
    seqres = []
    models = []            # list of (residues, atom_lines), one per MODEL
    cur_res, cur_lines = None, None
    with open(path) as f:
        for line in f:
            rec = line[:6].strip()
            if rec == "SEQRES":
                ch = line[11].strip()
                if ch == chain_want:
                    seqres.extend(r for r in line[19:70].split() if r in AA3)
            elif rec == "MODEL":
                cur_res, cur_lines = {}, []
            elif rec == "ENDMDL":
                if cur_res is not None:
                    models.append((cur_res, cur_lines))
                cur_res, cur_lines = None, None
            elif rec == "ATOM":
                ch = line[21].strip()
                if ch != chain_want:
                    continue
                try:
                    rs = int(line[22:26])
                except ValueError:
                    continue
                rn = line[17:20].strip()
                if rn not in AA3:
                    continue
                if cur_res is None:      # X-ray/cryo-EM: no MODEL records
                    cur_res, cur_lines = {}, []
                cur_lines.append(line)
                if line[12:16].strip() == "CA":
                    cur_res[rs] = (rn, (float(line[30:38]),
                                        float(line[38:46]),
                                        float(line[46:54])))
    if cur_res:
        models.append((cur_res, cur_lines))
    if not models:
        return "".join(AA3[r] for r in seqres), {}, []
    if len(models) == 1:
        residues, atom_lines = models[0]
    else:
        k = pick_medoid(models)
        residues, atom_lines = models[k]
        MODEL_PICK[path] = f"nmr_medoid: model {k + 1} of {len(models)}"
    seq = "".join(AA3[r] for r in seqres)
    return seq, residues, atom_lines


def nw_align(q, t, match=2, mismatch=-1, gap=-2):
    """Global alignment of query q onto target t. Returns list of
    (q_pos or None, t_pos or None) pairs (0-based positions)."""
    n, m = len(q), len(t)
    S = np.zeros((n + 1, m + 1), dtype=int)
    P = np.zeros((n + 1, m + 1), dtype=np.int8)  # 1=diag 2=up 3=left
    for i in range(1, n + 1):
        S[i, 0] = i * gap; P[i, 0] = 2
    for j in range(1, m + 1):
        S[0, j] = j * gap; P[0, j] = 3
    for i in range(1, n + 1):
        qi = q[i - 1]
        for j in range(1, m + 1):
            d = S[i - 1, j - 1] + (match if qi == t[j - 1] else mismatch)
            u = S[i - 1, j] + gap
            l = S[i, j - 1] + gap
            best = max(d, u, l)
            S[i, j] = best
            P[i, j] = 1 if best == d else (2 if best == u else 3)
    aln = []
    i, j = n, m
    while i > 0 or j > 0:
        p = P[i, j]
        if p == 1:
            aln.append((i - 1, j - 1)); i -= 1; j -= 1
        elif p == 2:
            aln.append((i - 1, None)); i -= 1
        else:
            aln.append((None, j - 1)); j -= 1
    return aln[::-1]


def write_renumbered(atom_lines, residues, resmap, out_path):
    """resmap: old resseq -> new resseq (input-sequence position, 1-based).
    Writes kept ATOM records with chain A and remapped numbering."""
    with open(out_path, "w") as f:
        for line in atom_lines:
            try:
                rs = int(line[22:26])
            except ValueError:
                continue
            if rs not in resmap:
                continue
            new = resmap[rs]
            f.write(line[:21] + "A" + f"{new:>4d}" + line[26:])
        f.write("END\n")


def kabsch_rmsd(a, b):
    ac = a - a.mean(axis=0)
    bc = b - b.mean(axis=0)
    cov = ac.T @ bc
    v, _, wt = np.linalg.svd(cov)
    d = np.sign(np.linalg.det(v @ wt))
    rot = v @ np.diag([1.0, 1.0, d]) @ wt
    return float(np.sqrt(((ac @ rot - bc) ** 2).sum(axis=1).mean()))


def tm_pair_score(a, b, n_iter=5):
    """TM-score between two equal-length arrays (same convention as eval)."""
    L = len(b)
    d0 = max(1.24 * (L - 15) ** (1.0 / 3.0) - 1.8, 0.5)
    aln = a - a.mean(axis=0)
    bb = b - b.mean(axis=0)
    v, _, wt = np.linalg.svd(aln.T @ bb)
    d = np.sign(np.linalg.det(v @ wt))
    rot = v @ np.diag([1.0, 1.0, d]) @ wt
    aln = aln @ rot + bb.mean(axis=0)
    for _ in range(n_iter):
        dd = np.sqrt(((aln - b) ** 2).sum(axis=1))
        w = 1.0 / (1.0 + (dd / d0) ** 2)
        w = w / w.sum()
        cm = (a * w[:, None]).sum(axis=0)
        cr = (b * w[:, None]).sum(axis=0)
        aa = a - cm
        bb = b - cr
        cov = (aa * w[:, None]).T @ bb
        v, _, wt = np.linalg.svd(cov)
        d = np.sign(np.linalg.det(v @ wt))
        rot = v @ np.diag([1.0, 1.0, d]) @ wt
        aln = aa @ rot + cr
    dd = np.sqrt(((aln - b) ** 2).sum(axis=1))
    return float((1.0 / (1.0 + (dd / d0) ** 2)).sum() / L)


def prepare_ref(pdb_id, chain, input_seq, cache_dir, out_path):
    """Download pdb, extract chain, align to input_seq, write renumbered PDB.
    Returns (n_matched, pos_map) where pos_map maps input position (1-based)
    -> CA xyz for positions covered by this reference."""
    os.makedirs(cache_dir, exist_ok=True)
    cache = os.path.join(cache_dir, f"{pdb_id.upper()}.pdb")
    if not os.path.exists(cache):
        fetch(f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb", dest=cache)
    seqres_seq, residues, atom_lines = parse_pdb_chain(cache, chain)
    if not residues:
        raise ValueError(f"{pdb_id} chain {chain}: no CA atoms")
    # chain sequence from ATOM records (sorted by resseq), aligned to input
    idx = sorted(residues)
    chain_seq = "".join(AA3[residues[i][0]] for i in idx)
    aln = nw_align(chain_seq, input_seq)
    resmap = {}
    pos_map = {}
    n_id = 0
    for q_pos, t_pos in aln:
        if q_pos is None:
            continue
        rs = idx[q_pos]
        if t_pos is not None:
            resmap[rs] = t_pos + 1
            pos_map[t_pos + 1] = residues[rs][1]
            if chain_seq[q_pos] == input_seq[t_pos]:
                n_id += 1
    write_renumbered(atom_lines, residues, resmap, out_path)
    ident = n_id / max(len(resmap), 1)
    return len(resmap), pos_map, seqres_seq, ident


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", required=True)
    ap.add_argument("--root", default=".", help="repo root (writes yamls/, panel_data/)")
    ap.add_argument("--only", default=None, help="comma-separated target names")
    ap.add_argument("--include-optional", action="store_true")
    args = ap.parse_args()

    only = set(args.only.split(",")) if args.only else None
    rows = list(csv.DictReader(open(args.targets)))
    # fetch_report.csv: one row per target, latest result (rewrite, not append)
    rep_path = os.path.join(args.root, "panel", "fetch_report.csv")
    rep_fields = ["name", "set", "L_input", "matched_ref1", "matched_ref2",
                  "shared", "id1", "id2", "ref_rmsd", "rmsd_table", "warning"]
    old = {}
    if os.path.exists(rep_path):
        for r in csv.DictReader(open(rep_path)):
            old[r["name"]] = r
    done_names = set()

    for row in rows:
        name = row["name"]
        if only and name not in only:
            continue
        if row["status"] != "ready":
            print(f"[skip] {name}: status={row['status']}")
            continue
        if row["optional"] == "1" and not args.include_optional:
            continue
        try:
            # ---- input sequence
            if row["uniprot"] and row.get("input_source", "uniprot") != "ref1_seqres":
                seq = uniprot_sequence(row["uniprot"])
            else:
                cache = os.path.join(args.root, "panel_data", "pdb_cache")
                os.makedirs(cache, exist_ok=True)
                p = os.path.join(cache, f"{row['ref1_pdb']}.pdb")
                if not os.path.exists(p):
                    fetch(f"https://files.rcsb.org/download/{row['ref1_pdb']}.pdb", dest=p)
                seq, res, _ = parse_pdb_chain(p, row["ref1_chain"])
                seq = seq  # SEQRES of ref1 chain
            if not seq:
                raise ValueError("empty input sequence")

            cache_dir = os.path.join(args.root, "panel_data", "pdb_cache")
            r1 = os.path.join(args.root, "panel_data", f"{name}_ref1.pdb")
            r2 = os.path.join(args.root, "panel_data", f"{name}_ref2.pdb")
            n1, m1, _, id1 = prepare_ref(row["ref1_pdb"], row["ref1_chain"], seq, cache_dir, r1)
            n2, m2, _, id2 = prepare_ref(row["ref2_pdb"], row["ref2_chain"], seq, cache_dir, r2)
            medoid_notes = [MODEL_PICK.get(os.path.join(cache_dir, f"{row[k].upper()}.pdb"), "")
                            for k in ("ref1_pdb", "ref2_pdb")]

            shared = sorted(set(m1) & set(m2))
            a = np.array([m1[k] for k in shared])
            b = np.array([m2[k] for k in shared])
            rr = kabsch_rmsd(a, b) if len(shared) >= 3 else float("nan")

            # ---- per-row validation gate (frozen, amendment 2026-10-05) ----
            warn = [n for n in medoid_notes if n]
            if min(id1, id2) < 0.95:
                warn.append(f"VALIDATION_FAIL: alignment identity {id1:.2f}/{id2:.2f} < 0.95")
            rt = row.get("rmsd_table", "")
            # Amendment 2026-10-05b (audit round 2, point 8): all-25-target check
            # shows AFsample2 table RMSD/TM values deviate from our metrics in
            # BOTH directions with no systematic offset — not a convertible
            # convention difference (their constructs/entries differ per row).
            # Table values are therefore INFORMATIONAL ONLY; the hard
            # transcription guard is the alignment-identity gate above.
            tm12 = tm_pair_score(a, b) if len(shared) >= 3 else float("nan")
            tm_t = row.get("tm_pair", "")
            if tm_t and tm12 == tm12 and abs(tm12 - float(tm_t)) > 0.10:
                warn.append(f"TM_INFO: ours {tm12:.2f} vs table {float(tm_t):.2f} (not a gate)")
            if len(shared) < 0.6 * len(seq):
                warn.append("low_shared_coverage")
            if rr == rr and rr < 3.5:
                warn.append("ref_rmsd<3.5: scaled basin threshold applies")
            if row["length"] and abs(len(seq) - int(row["length"])) > 5:
                warn.append(f"length mismatch (input {len(seq)} vs table {row['length']})")

            yd = os.path.join(args.root, "yamls")
            os.makedirs(yd, exist_ok=True)
            with open(os.path.join(yd, f"{name}.yaml"), "w") as f:
                f.write("version: 1\nsequences:\n  - protein:\n      id: A\n")
                f.write(f"      sequence: {seq}\n")

            old[name] = dict(zip(rep_fields, [name, row["set"], len(seq), n1, n2,
                                 len(shared), f"{id1:.3f}", f"{id2:.3f}",
                                 f"{rr:.2f}", rt, ";".join(warn)]))
            done_names.add(name)
            print(f"[ok] {name}: L={len(seq)} matched {n1}/{n2} shared={len(shared)} "
                  f"id={id1:.2f}/{id2:.2f} ref_rmsd={rr:.2f} "
                  f"{'WARN: ' + ';'.join(warn) if warn else ''}")
        except Exception as e:
            old[name] = dict(zip(rep_fields, [name, row["set"], "", "", "", "",
                                              "", "", "", "", f"ERROR: {e}"]))
            done_names.add(name)
            print(f"[FAIL] {name}: {e}", file=sys.stderr)

    with open(rep_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rep_fields)
        w.writeheader()
        for n, r in old.items():
            w.writerow(r)
    print(f"\nwrote {rep_path} ({len(old)} targets, latest row per target)")


if __name__ == "__main__":
    main()
