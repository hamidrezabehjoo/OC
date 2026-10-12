#!/usr/bin/env python3
"""Subsample a frozen panel MSA (.a3m) to a fixed depth, deterministically.

Used by the pre-registered MSA-depth-engineered boost arm (protocol
amendment 2026-10-05b): shallow MSAs collapse Boltz-2 diversity into the
dominant basin, giving the boost regime a guaranteed data point.

Rule (frozen): keep the query (first record) plus records at a constant
stride so the total depth equals --depth. No RNG.

    stride = max(1, (n_records - 1) // (depth - 1))
    kept   = query + records[0::stride], truncated to depth

Usage:
    python3 panel/msa_subsample.py in.a3m out.a3m --depth 32
"""

import argparse
import hashlib


def read_a3m(path):
    records = []
    header, seq = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq)))
                header, seq = line, []
            elif line.strip():
                seq.append(line.strip())
    if header is not None:
        records.append((header, "".join(seq)))
    return records


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("src")
    ap.add_argument("dst")
    ap.add_argument("--depth", type=int, default=32)
    args = ap.parse_args()

    records = read_a3m(args.src)
    if len(records) < 2:
        raise SystemExit(f"{args.src}: only {len(records)} record(s), nothing to do")
    query, rest = records[0], records[1:]
    if args.depth < 2:
        raise SystemExit("--depth must be >= 2 (query + >=1 homolog)")
    stride = max(1, len(rest) // (args.depth - 1))
    kept = [query] + rest[::stride][: args.depth - 1]
    with open(args.dst, "w") as f:
        for h, s in kept:
            f.write(h + "\n")
            f.write(s + "\n")
    md5 = hashlib.md5(open(args.dst, "rb").read()).hexdigest()[:12]
    print(f"{args.src}: {len(records)} records -> {args.dst}: {len(kept)} "
          f"(stride {stride}) md5 {md5}")


if __name__ == "__main__":
    main()
