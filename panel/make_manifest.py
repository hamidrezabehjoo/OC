#!/usr/bin/env python3
"""Write or verify the panel freeze manifest (amendment 2026-10-05c).

    python3 panel/make_manifest.py            # write panel/MANIFEST.md5
    python3 panel/make_manifest.py --verify   # check current files against it

The manifest pins md5s of every input the screen depends on: targets.csv,
the panel scripts, the frozen protocol, prepared references
(panel_data/*_ref*.pdb), boltz yamls, frozen MSAs (panel_msa/*.a3m,
when present), and the patched controller file diffusionv2.py (the
NU_ALPHA hook + numerics guards — the method under test, pinned
2026-10-07 after the freeze-coverage gap was noticed). On the GPU
machine, run with --verify instead of re-fetching: a re-fetch hits
RCSB again and could change references.
"""

import argparse
import glob
import hashlib
import os
import sys

MANIFEST = os.path.join("panel", "MANIFEST.md5")

PATTERNS = [
    "panel/targets.csv",
    "panel/*.py",
    "panel/*.sh",
    "panel/prospective_test_panel_protocol.md",
    "panel_data/*_ref1.pdb",
    "panel_data/*_ref2.pdb",
    "yamls/*.yaml",
    "panel_msa/*.a3m",
    "diffusionv2.py",
]


def md5(path):
    return hashlib.md5(open(path, "rb").read()).hexdigest()


def collect():
    files = []
    for pat in PATTERNS:
        files.extend(glob.glob(pat))
    return sorted(set(files))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--verify", action="store_true")
    args = ap.parse_args()

    if not args.verify:
        with open(MANIFEST, "w") as f:
            for p in collect():
                f.write(f"{md5(p)}  {p}\n")
        print(f"wrote {MANIFEST} ({len(collect())} files)")
        return

    if not os.path.exists(MANIFEST):
        sys.exit(f"{MANIFEST} not found")
    bad = missing = 0
    for line in open(MANIFEST):
        line = line.rstrip("\n")
        if not line:
            continue
        want, path = line.split(None, 1)
        path = path.strip()
        if not os.path.exists(path):
            print(f"MISSING  {path}")
            missing += 1
        elif md5(path) != want:
            print(f"CHANGED  {path}")
            bad += 1
    if bad or missing:
        sys.exit(f"VERIFY FAILED: {bad} changed, {missing} missing")
    print("manifest OK: all pinned files match")


if __name__ == "__main__":
    main()
