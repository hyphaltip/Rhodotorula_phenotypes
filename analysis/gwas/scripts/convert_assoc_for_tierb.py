#!/usr/bin/env python3
"""Convert GEMMA's raw per-trait .assoc.txt (Tier A rebuild output) into the
`{prefix}_{trait}_assoc.csv.gz` format tierb_set_tests.py expects (rs, chr, ps,
af, beta, se, p_wald columns; matches the naming/shape of the original run's
results/gwas/tierA_summary/{gwas,gwasc}_<trait>_assoc.csv.gz files).
"""
import argparse
import glob
import os
import re

import pandas as pd

COLS = ["rs", "chr", "ps", "af", "beta", "se", "p_wald"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--assoc-dir", required=True, help="dir with <panel>_<trait>.assoc.txt")
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    files = sorted(glob.glob(os.path.join(args.assoc_dir, "*.assoc.txt")))
    assert files, f"no .assoc.txt files found in {args.assoc_dir}"

    for f in files:
        stem = re.sub(r"\.assoc\.txt$", "", os.path.basename(f))
        df = pd.read_csv(f, sep="\t", usecols=COLS)
        out = os.path.join(args.out_dir, f"{stem}_assoc.csv.gz")
        df.to_csv(out, index=False, compression="gzip")
        print(f"  {stem}: {len(df)} SNPs -> {out}")

    print(f"Converted {len(files)} assoc files to {args.out_dir}")


if __name__ == "__main__":
    main()
