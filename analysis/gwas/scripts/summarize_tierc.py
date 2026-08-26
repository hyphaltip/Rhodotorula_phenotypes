#!/usr/bin/env python3
"""Summarize Tier C BSLMM output: PVE/PGE/n_gamma posterior summaries (from
.hyp.txt.gz) and per-SNP posterior inclusion probability (PIP, from
.gamma.txt.gz) -- matches the format of
analysis/ideas/2026-08-15-color-phenotype-space/results/gwas/tierC_summary/tierc_bslmm_summary.csv.

GEMMA's .gamma.txt.gz has one row per recorded MCMC sample, columns s0..s(W-1)
(W = the max n_gamma seen across all samples), holding the 0-indexed SNP
positions (row numbers into the .bim) included in that sample's sparse
component. Only the first n_gamma[row] entries per row are real -- later
columns are zero-padding, and a genuine SNP index of 0 is indistinguishable
from padding by value alone, so n_gamma (from the paired .hyp.txt.gz row) MUST
be used to truncate each row before counting. PIP for a SNP = (# samples where
its index appears in the first n_gamma entries) / (total samples).
"""
import argparse
import gzip

import numpy as np
import pandas as pd


def load_bim_rsids(bim_path: str) -> list[str]:
    ids = []
    with open(bim_path) as f:
        for line in f:
            parts = line.split()
            ids.append(parts[1])
    return ids


def compute_pip(hyp_path: str, gamma_path: str, bim_rsids: list[str]) -> pd.DataFrame:
    hyp = pd.read_csv(hyp_path, sep="\t", index_col=False)
    hyp.columns = [c.strip() for c in hyp.columns]
    n_gamma = hyp["n_gamma"].to_numpy(dtype=int)

    counts: dict[int, int] = {}
    n_samples = 0
    with gzip.open(gamma_path, "rt") as f:
        header = f.readline()  # s0 s1 ...
        for i, line in enumerate(f):
            row = line.split()
            if not row:
                continue
            k = int(n_gamma[i]) if i < len(n_gamma) else 0
            for tok in row[:k]:
                idx = int(tok)
                counts[idx] = counts.get(idx, 0) + 1
            n_samples += 1

    rows = []
    for idx, cnt in counts.items():
        pip = cnt / n_samples if n_samples else 0.0
        rs = bim_rsids[idx] if 0 <= idx < len(bim_rsids) else f"UNKNOWN_IDX_{idx}"
        rows.append({"bim_index": idx, "rs": rs, "pip": pip, "n_samples_included": cnt})
    df = pd.DataFrame(rows).sort_values("pip", ascending=False).reset_index(drop=True)
    df.attrs["n_mcmc_samples"] = n_samples
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tierc-dir", required=True)
    ap.add_argument("--bim", required=True, help=".bim of the bfile BSLMM was run on (index-matched)")
    ap.add_argument("--traits", nargs="+", required=True)
    ap.add_argument("--out-summary", required=True)
    ap.add_argument("--out-pip-dir", required=True)
    ap.add_argument("--pip-threshold", type=float, default=0.5,
                     help="report loci at or above this PIP in the summary's top-loci column")
    args = ap.parse_args()

    import os
    os.makedirs(args.out_pip_dir, exist_ok=True)
    os.makedirs(os.path.dirname(args.out_summary), exist_ok=True)

    bim_rsids = load_bim_rsids(args.bim)

    rows = []
    for trait in args.traits:
        hyp_path = os.path.join(args.tierc_dir, f"bslmm_{trait}.hyp.txt.gz")
        gamma_path = os.path.join(args.tierc_dir, f"bslmm_{trait}.gamma.txt.gz")
        hyp = pd.read_csv(hyp_path, sep="\t", index_col=False)
        hyp.columns = [c.strip() for c in hyp.columns]

        pve = hyp["pve"].to_numpy()
        pge = hyp["pge"].to_numpy()
        ngam = hyp["n_gamma"].to_numpy()

        pip_df = compute_pip(hyp_path, gamma_path, bim_rsids)
        pip_df.to_csv(os.path.join(args.out_pip_dir, f"{trait}_pip.csv"), index=False)

        top_loci = pip_df[pip_df["pip"] >= args.pip_threshold]
        top_desc = "; ".join(f"{r.rs} PIP={r.pip:.2f}" for r in top_loci.itertuples()) or "none >= threshold"
        n_clusters_pip_gt = len(top_loci)

        rows.append({
            "trait": trait,
            "n_mcmc_samples": pip_df.attrs["n_mcmc_samples"],
            "PVE_med": round(float(np.median(pve)), 3),
            "PVE_lo95": round(float(np.percentile(pve, 2.5)), 3),
            "PVE_hi95": round(float(np.percentile(pve, 97.5)), 3),
            "PGE_med": round(float(np.median(pge)), 3),
            "n_gamma_med": round(float(np.median(ngam)), 1),
            f"n_loci_PIP_ge_{args.pip_threshold}": n_clusters_pip_gt,
            "top_PIP_loci": top_desc,
        })
        print(f"{trait}: PVE={np.median(pve):.3f} [{np.percentile(pve,2.5):.3f},{np.percentile(pve,97.5):.3f}]  "
              f"PGE={np.median(pge):.3f}  n_gamma_med={np.median(ngam):.1f}  "
              f"n_loci(PIP>={args.pip_threshold})={n_clusters_pip_gt}")

    out = pd.DataFrame(rows)
    out.to_csv(args.out_summary, index=False)
    print(f"\nWrote {args.out_summary}")
    print(f"Wrote per-trait PIP tables to {args.out_pip_dir}/")


if __name__ == "__main__":
    main()
