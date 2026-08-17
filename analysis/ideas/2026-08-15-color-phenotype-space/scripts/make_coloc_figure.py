#!/usr/bin/env python3
"""dxy/Fst co-localization figure:
  [A] per-window mean dxy vs mean Fst, highlighting windows that contain a
      Tier-A FDR SNP (any trait); anchor loci annotated
  [B] bar chart: observed vs expected (genome-wide rate) fraction of GWAS-windows
      that are high-dxy / high-Fst
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
    ap.add_argument("--coloc-final", required=True)
    ap.add_argument("--anchors", required=True)
    ap.add_argument("--dxy", required=True)
    ap.add_argument("--fst", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()
    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)
    f = pd.read_csv(args.coloc_final)
    a = pd.read_csv(args.anchors)
    dxy = pd.read_csv(args.dxy, sep="\t")
    fst = pd.read_csv(args.fst, sep="\t")
    win = pd.DataFrame({
        "chromosome": dxy.groupby("chromosome")["avg_dxy"].count().index,
    })
    win = win.merge(dxy.groupby(["chromosome", "window_pos_1"])["avg_dxy"].mean()
                    .reset_index(), on="chromosome", how="inner")
    win = win.merge(fst.groupby(["chromosome", "window_pos_1"])["avg_hudson_fst"]
                    .mean().reset_index(), on=["chromosome", "window_pos_1"])
    win = win.rename(columns={"avg_dxy": "mean_dxy", "avg_hudson_fst": "mean_fst"})
    win["high_dxy"] = win["mean_dxy"] >= win["mean_dxy"].quantile(0.80)
    win["high_fst"] = win["mean_fst"] >= win["mean_fst"].quantile(0.80)

    # mark windows that contain >=1 FDR SNP (f has chr_s + ps, win size 100kb)
    def in_win(row):
        m = win[(win["chromosome"] == row["chr_s"]) &
                (win["window_pos_1"] <= row["ps"]) &
                (row["ps"] < win["window_pos_1"] + 100000)]
        return m.index[0] if len(m) else -1
    hitidx = set(int(x) for x in f.apply(in_win, axis=1) if x != -1)
    win["has_gwas"] = [i in hitidx for i in range(len(win))]

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    ax = axes[0]
    o = ~win["has_gwas"]
    ax.scatter(win.loc[o, "mean_fst"], win.loc[o, "mean_dxy"], s=12, alpha=.6,
               c="#9ecae1", edgecolors="none", label="windows w/o FDR SNP")
    g = win.columns
    ax.scatter(win.loc[~o, "mean_fst"], win.loc[~o, "mean_dxy"], s=12, c="#3182bd",
               edgecolors="none", label="Windows w/ FDR SNP")
    ax.scatter(a["mean_fst"], a["mean_dxy"], s=70, facecolors="none",
               edgecolors="#e31a1c", linewidths=1.4, label="GWAS anchor loci")
    for _, r in a.iterrows():
        ax.annotate(r["trait"], (r["mean_fst"], r["mean_dxy"]),
                    xytext=(5, 5), textcoords="offset points", fontsize=7,
                    color="#e31a1c")
    ax.set_xlabel("mean Fst (15 pop pairs)"); ax.set_ylabel("mean dxy (15 pop pairs)")
    ax.set_title("(A) window divergence vs GWAS co-localization")
    ax.legend(fontsize=7)

    ax = axes[1]
    obs = win.loc[win["has_gwas"], ["high_dxy", "high_fst"]].mean()
    exp = win[["high_dxy", "high_fst"]].mean()
    x = np.arange(2)
    ax.bar(x - 0.18, exp, width=0.36, label="genome-wide rate", color="#969696")
    ax.bar(x + 0.18, obs, width=0.36, label="GWAS windows", color="#3182bd")
    ax.set_xticks(x); ax.set_xticklabels(["high-dxy", "high-Fst"])
    ax.set_ylabel("fraction of windows")
    ax.set_title("(B) GWAS windows vs genome rate")
    ax.legend(fontsize=7)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(out / f"coloc_dxy_fst.{ext}", dpi=150)
    print("wrote coloc_dxy_fst.{png,pdf}")

if __name__ == "__main__":
    main()
