#!/usr/bin/env python3
"""Tier B set-test figure:
  [A] QQ-style calibration for burden_p & skat_p (per-set overlays)
  [B] genomic windows sorted by skat_p (with burden & min_p), marking locus of
      any window with skat_p<0.01 and FDR thresholds
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
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()
    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8))
    markers = ["o", "s"]
    for i, f in enumerate(args.tierb):
        pref = Path(f).stem.replace("tierb_settests_", "")
        df = pd.read_csv(f)
        for col, style, ax in [("burden_p", "-", axes[0]), ("skat_p", "-", axes[0])]:
            p = df[col].dropna().to_numpy()
            p = p[np.isfinite(p) & (p > 0)]
            p = np.sort(p)
            n = len(p)
            exp = (np.arange(1, n + 1) - 0.5) / n
            ax.plot(-np.log10(exp), -np.log10(p),
                    markers[i] + style, ms=3, lw=0.6, alpha=.8,
                    label=f"{pref} {col}")
    axes[0].plot([0, 5], [0, 5], "k--", lw=0.8)
    axes[0].set_xlabel("expected -log10 p"); axes[0].set_ylabel("observed -log10 p")
    axes[0].set_title("(A) set-test calibration")
    axes[0].legend(fontsize=7, ncol=2)
    axes[0].set_xlim(0, 4.5); axes[0].set_ylim(0, 4.5)

    # B: per-set skat window ranking
    ax = axes[1]
    for i, f in enumerate(args.tierb):
        pref = Path(f).stem.replace("tierb_settests_", "")
        df = pd.read_csv(f)
        d = df.sort_values("skat_p").reset_index(drop=True)
        ax.plot(d.index + 1, -np.log10(d["skat_p"]), markers[i], ms=2.5,
                alpha=.8, label=f"{pref} (skat)")
        # annotate best
        q = d.iloc[0]
        ax.annotate(f"{q['scaffold']}:{q['win_start']}", (1, -np.log10(q['skat_p'])),
                    xytext=(8, 8), textcoords="offset points", fontsize=7)
    ax.axhline(-np.log10(0.05), color="r", ls="--", lw=0.8)
    ax.set_xlabel("window rank (by skat_p)"); ax.set_ylabel("-log10 p")
    ax.set_yscale("log")
    ax.set_title("(B) SKAT across windows")
    ax.legend(fontsize=7)
    ax.text(.98, .98, "red = p 0.05", transform=ax.transAxes, ha="right", va="top", fontsize=7)

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
