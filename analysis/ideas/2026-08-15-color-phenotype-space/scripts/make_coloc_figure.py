#!/usr/bin/env python3
"""dxy/Fst co-localization figure:
  [A] per-window mean dxy vs mean Fst, highlighting windows that contain a
      Tier-A FDR SNP (any trait); 80th-pct high-dxy / high-Fst threshold lines;
      anchor loci annotated
  [B] observed vs expected (genome-wide rate) fraction of GWAS-windows that are
      high-dxy / high-Fst, with Wilson 95% CIs and Fisher-exact p-values
Saves PNG + PDF.
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact

def wilson_ci(k, n, z=1.96):
    """Wilson score 95% CI for a binomial proportion."""
    if n == 0:
        return 0.0, 0.0, 0.0
    p_hat = k / n
    denom = 1 + z**2 / n
    centre = (p_hat + z**2 / (2 * n)) / denom
    half = z * np.sqrt(p_hat * (1 - p_hat) / n + z**2 / (4 * n**2)) / denom
    return p_hat, max(0.0, centre - half), min(1.0, centre + half)

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
    # thresholds MATCH the enrichment computation (80th percentile across windows)
    thr_dxy = win["mean_dxy"].quantile(0.80)
    thr_fst = win["mean_fst"].quantile(0.80)
    win["high_dxy"] = win["mean_dxy"] >= thr_dxy
    win["high_fst"] = win["mean_fst"] >= thr_fst

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
               edgecolors="none", label="windows w/ FDR SNP")
    ax.axvline(thr_fst, color="k", ls=":", lw=0.9)
    ax.axhline(thr_dxy, color="k", ls=":", lw=0.9)
    ax.text(0.98, .98, f"dotted = 80th pct\n(high-dxy {thr_dxy:.3f}, high-Fst {thr_fst:.3f})",
            transform=ax.transAxes, ha="right", va="top", fontsize=7,
            bbox=dict(fc="w", ec="0.8", lw=0.6))
    ax.scatter(a["mean_fst"], a["mean_dxy"], s=70, facecolors="none",
               edgecolors="#e31a1c", linewidths=1.4, label="GWAS anchor loci")
    # de-collide anchor labels sharing the same window (chroma & AUC_10 both scaffold_10)
    used = {}
    for _, r in a.sort_values("trait").iterrows():
        x0 = float(r["mean_fst"]); y0 = float(r["mean_dxy"])
        key = (round(x0, 5), round(y0, 5))
        if key in used:                 # same window -> stack labels
            used[key] += 1
            y0 = y0 - used[key] * 0.006
        else:
            used[key] = 0
        ax.annotate(r["trait"], (r["mean_fst"], r["mean_dxy"]),
                    xytext=(5, 5 - used[key] * 8), textcoords="offset points",
                    fontsize=7, color="#e31a1c")
    ax.set_xlabel("mean Fst (15 pop pairs)"); ax.set_ylabel("mean dxy (15 pop pairs)")
    ax.set_title("(A) window divergence vs GWAS co-localization")
    ax.legend(fontsize=7, loc="lower right")

    ax = axes[1]
    n_gwas = int(win["has_gwas"].sum())
    xt = [("high-dxy", "high_dxy"), ("high-Fst", "high_fst")]
    for i, (name, col) in enumerate(xt):
        # genome-wide rate across ALL windows
        k_all = int(win[col].sum()); n_all = len(win)
        ephat, elo, ehi = wilson_ci(k_all, n_all)
        # GWAS windows
        k_gw = int((win["has_gwas"] & win[col]).sum()); n_gw = n_gwas
        phat, lo, hi = wilson_ci(k_gw, n_gw)
        # Fisher: gwas high/low vs not-gwas high/low
        g_high, g_low = k_gw, n_gw - k_gw
        ng_high, ng_low = k_all - k_gw, (n_all - n_gw) - (k_all - k_gw)
        table = [[g_high, g_low], [ng_high, ng_low]]
        orr, pv = fisher_exact(table)   # alternative='two-sided'
        x = np.arange(len(xt))
        ax.bar(x[i] - 0.18, ephat, width=0.36, color="#969696",
               yerr=[[ephat - elo], [ehi - ephat]], capsize=3,
               label="genome-wide rate" if i == 0 else None,
               error_kw=dict(lw=0.9))
        ax.bar(x[i] + 0.18, phat, width=0.36, color="#3182bd",
               yerr=[[phat - lo], [hi - phat]], capsize=3,
               label="GWAS windows" if i == 0 else None,
               error_kw=dict(lw=0.9))
        ax.text(x[i] + 0.18, hi + 0.02, f"n={k_gw}/{n_gw}", ha="center", fontsize=7)
        ax.text(x[i] + 0.18, phat / 2 if phat else 0.05,
                f"OR={orr:.2f}\np={pv:.3f}", ha="center", va="center",
                fontsize=7, color="#0b3d91")
    ax.set_xticks(np.arange(len(xt))); ax.set_xticklabels(xt)
    ax.set_ylabel("fraction of windows")
    ax.set_ylim(0, 0.6)
    ax.set_title("(B) GWAS windows vs genome rate (95% Wilson CI)")
    ax.legend(fontsize=7)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(out / f"coloc_dxy_fst.{ext}", dpi=150)
    print("wrote coloc_dxy_fst.{png,pdf}")

if __name__ == "__main__":
    main()
