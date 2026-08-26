#!/usr/bin/env python3
"""Summarize Tier A GEMMA scans: genomic inflation (lambda), top hit, and
BH-FDR(q=0.05) significant-SNP count per trait -- matches the summary format in
analysis/ideas/2026-08-15-color-phenotype-space/PROGRESS.md sections 6 and 9.

Also reports a Meff (effective number of independent tests) PROXY and the
corresponding Bonferroni threshold, per the 2026-08-25 quant-genetics consult
(D-16): BH-FDR alone doesn't flag that many "significant" SNPs in a near-clonal
panel are near-duplicate observations from one LD block. A full LD-eigenvalue
Meff (Li & Ji / Galwey) over ~500k unpruned SNPs is not computed here (too
costly); instead --n-pruned-snps (the LD-pruned marker count used for this
panel's kinship, e.g. 29,453 for gwas/213) is reported as a standard rough
Meff proxy -- pruning already targeted r^2<0.2 approximate independence, so its
SNP count approximates "how many independent tests" exist genome-wide. This is
a documented approximation, not an exact Meff -- treat n_fdr05 as the primary
significance call and meff_proxy/bonferroni_meff as sensitivity context.
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
    ap.add_argument("--prefix", default="gwas", help="file prefix, e.g. 'gwas' or 'gwasc'")
    ap.add_argument("--out-summary", required=True)
    ap.add_argument("--out-fdr-dir", required=True)
    ap.add_argument("--n-pruned-snps", type=int, default=None,
                     help="LD-pruned marker count for this panel's kinship -- used as a Meff proxy")
    args = ap.parse_args()

    os.makedirs(args.out_fdr_dir, exist_ok=True)
    os.makedirs(os.path.dirname(args.out_summary), exist_ok=True)

    files = sorted(glob.glob(os.path.join(args.assoc_dir, f"{args.prefix}_*.assoc.txt")))
    assert files, f"no assoc files found in {args.assoc_dir} matching {args.prefix}_*.assoc.txt"

    bonf_meff = (0.05 / args.n_pruned_snps) if args.n_pruned_snps else None

    rows = []
    for f in files:
        trait = re.sub(rf"^{args.prefix}_|\.assoc\.txt$", "", os.path.basename(f))
        df = pd.read_csv(f, sep="\t")
        df = df[df["p_wald"].notna() & (df["p_wald"] > 0)]
        n = len(df)
        lam = lambda_gc(df["p_wald"].to_numpy())
        top = df.loc[df["p_wald"].idxmin()]
        sig_mask = bh_fdr(df["p_wald"].to_numpy(), q=0.05)
        n_fdr = int(sig_mask.sum())
        n_bonf_meff = int((df["p_wald"] < bonf_meff).sum()) if bonf_meff else None
        df.loc[sig_mask].sort_values("p_wald").to_csv(
            os.path.join(args.out_fdr_dir, f"{args.prefix}_{trait}_fdr05.csv"), index=False
        )
        rows.append({
            "trait": trait,
            "n_snps": n,
            "lambda_gc": round(lam, 3),
            "n_fdr05": n_fdr,
            "meff_proxy": args.n_pruned_snps,
            "bonferroni_meff_thresh": bonf_meff,
            "n_sig_bonferroni_meff": n_bonf_meff,
            "top_snp": top["rs"],
            "top_chr": top["chr"],
            "top_ps": top["ps"],
            "top_p_wald": top["p_wald"],
            "top_beta": top["beta"],
        })
        extra = f" n_sig_bonf(meff={args.n_pruned_snps})={n_bonf_meff}" if bonf_meff else ""
        print(f"{trait}: n={n} lambda={lam:.3f} n_FDR05={n_fdr}{extra} top={top['rs']} p={top['p_wald']:.3e}")

    out = pd.DataFrame(rows).sort_values("trait")
    out.to_csv(args.out_summary, index=False)
    print(f"\nWrote {args.out_summary}")
    print(f"Wrote per-trait FDR05 SNP lists to {args.out_fdr_dir}/")


if __name__ == "__main__":
    main()
