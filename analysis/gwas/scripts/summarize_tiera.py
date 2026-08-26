#!/usr/bin/env python3
"""Summarize Tier A GEMMA scans: genomic inflation (lambda), top hit, and
BH-FDR(q=0.05) significant-SNP count per trait -- matches the summary format in
analysis/ideas/2026-08-15-color-phenotype-space/PROGRESS.md sections 6 and 9.
"""
import argparse
import glob
import os
import re

import numpy as np
import pandas as pd
from scipy.stats import chi2


def lambda_gc(pvals: np.ndarray) -> float:
    chisq = chi2.isf(pvals, df=1)
    return float(np.median(chisq) / chi2.median(df=1))


def bh_fdr(pvals: np.ndarray, q: float = 0.05) -> np.ndarray:
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order]
    thresh = (np.arange(1, n + 1) / n) * q
    below = ranked <= thresh
    if not below.any():
        return np.zeros(n, dtype=bool)
    max_i = np.max(np.where(below)[0])
    cutoff_p = ranked[max_i]
    return pvals <= cutoff_p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--assoc-dir", required=True)
    ap.add_argument("--out-summary", required=True)
    ap.add_argument("--out-fdr-dir", required=True)
    args = ap.parse_args()

    os.makedirs(args.out_fdr_dir, exist_ok=True)
    os.makedirs(os.path.dirname(args.out_summary), exist_ok=True)

    files = sorted(glob.glob(os.path.join(args.assoc_dir, "gwas_*.assoc.txt")))
    assert files, f"no assoc files found in {args.assoc_dir}"

    rows = []
    for f in files:
        trait = re.sub(r"^gwas_|\.assoc\.txt$", "", os.path.basename(f))
        df = pd.read_csv(f, sep="\t")
        df = df[df["p_wald"].notna() & (df["p_wald"] > 0)]
        n = len(df)
        lam = lambda_gc(df["p_wald"].to_numpy())
        top = df.loc[df["p_wald"].idxmin()]
        sig_mask = bh_fdr(df["p_wald"].to_numpy(), q=0.05)
        n_fdr = int(sig_mask.sum())
        df.loc[sig_mask].sort_values("p_wald").to_csv(
            os.path.join(args.out_fdr_dir, f"{trait}_fdr05.csv"), index=False
        )
        rows.append({
            "trait": trait,
            "n_snps": n,
            "lambda_gc": round(lam, 3),
            "n_fdr05": n_fdr,
            "top_snp": top["rs"],
            "top_chr": top["chr"],
            "top_ps": top["ps"],
            "top_p_wald": top["p_wald"],
            "top_beta": top["beta"],
        })
        print(f"{trait}: n={n} lambda={lam:.3f} n_FDR05={n_fdr} top={top['rs']} p={top['p_wald']:.3e}")

    out = pd.DataFrame(rows).sort_values("trait")
    out.to_csv(args.out_summary, index=False)
    print(f"\nWrote {args.out_summary}")
    print(f"Wrote per-trait FDR05 SNP lists to {args.out_fdr_dir}/")


if __name__ == "__main__":
    main()
