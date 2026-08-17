#!/usr/bin/env python3
"""Tier B set-test figure:
  [A] QQ-style calibration for burden_p & skat_p (per-set overlays)
  [B] genomic windows sorted by skat_p (with burden & min_p), marking locus of
      any window with skat_p<0.01 and FDR thresholds; exact-MC verified
      top-50 skat p overlaid if provided
  [C] FDR hit counts bar for burden/skat/min_p (highdxy subset)
Saves PNG + PDF.
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tierb", nargs="+", required=True)
    ap.add_argument("--mcver", nargs="+", default=None,
                    help="tierb_skat_mcver_{gwas,gwasc}.csv (optional, overlaid on B)")
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()
    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)

    mc_df = {}
    if args.mcver:
        for f in args.mcver:
            pref = Path(f).stem.replace("tierb_skat_mcver_", "")
            mc_df[pref] = pd.read_csv(f)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8))
    markers = ["o", "s"]
    for i, f in enumerate(args.tierb):
        pref = Path(f).stem.replace("tierb_settests_", "")
        df = pd.read_csv(f)
        for col in ["burden_p", "skat_p"]:
            p = df[col].dropna().to_numpy()
            p = p[np.isfinite(p) & (p > 0)]
            p = np.sort(p)
            n = len(p)
            exp = (np.arange(1, n + 1) - 0.5) / n
            axes[0].plot(-np.log10(exp), -np.log10(p),
                         markers[i] + "-", ms=3, lw=0.6, alpha=.8,
                         label=f"{pref} {col}")
    axes[0].plot([0, 5], [0, 5], "k--", lw=0.8)
    axes[0].set_xlabel("expected -log10 p"); axes[0].set_ylabel("observed -log10 p")
    axes[0].set_title("(A) set-test calibration")
    axes[0].legend(fontsize=7, ncol=2)
    axes[0].set_xlim(0, 4.5); axes[0].set_ylim(0, 4.5)

    # B: per-set skat window ranking (linear -log10 scale)
    ax = axes[1]
    for i, f in enumerate(args.tierb):
        pref = Path(f).stem.replace("tierb_settests_", "")
        df = pd.read_csv(f)
        d = df.sort_values("skat_p").reset_index(drop=True)
        yrank = np.where(d["skat_p"] > 0, -np.log10(d["skat_p"].clip(lower=1e-300)), 8.0)
        ax.plot(d.index + 1, yrank, markers[i], ms=2.5, alpha=.8,
                label=f"{pref} (SKAT: moment-approx)")
        q = d.iloc[0]
        ax.annotate(f"{q['scaffold']}:{q['win_start']}",
                    (1, float(yrank[0])),
                    xytext=(8, 8), textcoords="offset points", fontsize=7)
        if pref in mc_df:
            m = mc_df[pref]
            rank_of = d[["trait", "scaffold", "win_start"]].reset_index()
            j = rank_of.merge(m[["trait", "scaffold", "win_start", "skat_p_mcver"]],
                              on=["trait", "scaffold", "win_start"], how="inner")
            if len(j):
                jmc = j.sort_values("index")
                mxp = np.where(jmc["skat_p_mcver"] > 0,
                               -np.log10(jmc["skat_p_mcver"].clip(lower=1e-300)), 8.0)
                ax.plot(jmc["index"] + 1, mxp, "x", ms=4,
                        color="#7a0177", alpha=.9,
                        label=f"{pref} (exact MC, top-50)")
    ax.axhline(-np.log10(0.05), color="r", ls="--", lw=0.8)
    ax.axhline(-np.log10(0.001), color="r", ls=":", lw=0.8)
    ax.set_xlabel("window rank (by skat_p)"); ax.set_ylabel("-log10 p")
    ax.set_title("(B) SKAT across windows")
    ax.legend(fontsize=7)
    ax.text(.98, .98, "dashed = p 0.05; dotted = p 0.001",
            transform=ax.transAxes, ha="right", va="top", fontsize=7)

    # C: FDR highdxy hit counts
    ax = axes[2]
    cols = ["fdr_burden_p_highdxy", "fdr_skat_p_highdxy", "fdr_min_p_highdxy"]
    labels = ["burden", "SKAT", "min-p"]
    for i, f in enumerate(args.tierb):
        pref = Path(f).stem.replace("tierb_settests_", "")
        df = pd.read_csv(f)
        cnts = [int((df[c].dropna() < 0.05).sum()) for c in cols]
        x = np.arange(len(cols)) + i * 0.4
        ax.bar(x, cnts, width=0.4, label=pref)
        for xi, c in zip(x, cnts):
            ax.text(xi, c + 1, str(c), ha="center", fontsize=9)
    ax.set_xticks(np.arange(len(cols)) + 0.2); ax.set_xticklabels(labels)
    ax.set_ylabel("FDR(q<0.05) sig windows (high-dxy subset)")
    ax.set_title("(C) FDR hits in high-dxy windows")
    ax.legend(fontsize=7)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(out / f"tierB_settests.{ext}", dpi=150)
    print("wrote tierB_settests.{png,pdf}")

if __name__ == "__main__":
    main()
