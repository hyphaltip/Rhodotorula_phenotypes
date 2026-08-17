#!/usr/bin/env python3
"""Manhattan + QQ plots for next-gen Tier-A GWAS results.

Reads {gwas,gwasc}_<trait>_assoc.csv.gz from results/gwas/tierA_summary/ and
multi_{gwas,gwasc}_{color,copper}_fisher_top500.csv, produces per-trait
manhattan + QQ PNGs plus a combined panel for the key traits.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
import os

HERE = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
SRC = os.path.join(HERE, "results/gwas/tierA_summary")
OUT = os.path.join(HERE, "results/gwas/figures")
os.makedirs(OUT, exist_ok=True)

KEY_TRAITS = ["chroma", "sat", "bright", "AUC_10", "AUC_20", "AUC_30",
              "resilience_30", "cu_dose_slope"]
SETS = ["gwas", "gwasc"]
GW = 5e-8
SUGG = 5e-6


def manhattan(ax, df, title, pcol="p_wald"):
    """df: chr(int), ps, pcol. x = genomic position along chroms."""
    chroms = sorted(df["chr"].unique())
    offsets = {}
    cum = 0
    for c in chroms:
        offsets[c] = cum
        cum += df.loc[df["chr"] == c, "ps"].max() + 1
    xs = df["chr"].map(offsets) + df["ps"]
    colors = ["#4477AA", "#EE6677"]
    for i, c in enumerate(chroms):
        m = df["chr"] == c
        ax.scatter(xs[m], -np.log10(df.loc[m, pcol]), s=6,
                   c=colors[i % 2], alpha=0.6, linewidths=0)
    ax.axhline(-np.log10(GW), color="red", lw=1, ls="--", label="5e-8")
    ax.axhline(-np.log10(SUGG), color="orange", lw=1, ls=":", label="5e-6")
    ax.set_xlabel("Scaffold (genomic position)")
    ax.set_ylabel(r"$-\log_{10}$ p")
    ax.set_title(title, fontsize=9)
    # chrom labels at midpoints
    mids = []
    lbls = []
    prev_end = 0
    for c in chroms:
        end = offsets[c] + df.loc[df["chr"] == c, "ps"].max() + 1
        mids.append((prev_end + end) / 2)
        lbls.append(str(c))
        prev_end = end
    ax.set_xticks(mids)
    ax.set_xticklabels(lbls, fontsize=6)
    ax.legend(fontsize=6, loc="upper right")


def qq(ax, df, title):
    p = df["p_wald"].dropna().values
    p = p[(p > 0) & (p < 1)]
    n = len(p)
    exp = -np.log10(np.linspace(1 / (n + 1), n / (n + 1), n))
    obs = -np.log10(np.sort(p))
    ax.scatter(exp, obs, s=4, alpha=0.5, c="#4477AA", linewidths=0)
    lo, hi = max(exp), max(obs)
    ax.plot([0, max(lo, hi)], [0, max(lo, hi)], "k--", lw=1)
    # lambda
    chi2 = stats.chi2.isf(p, df=1)
    lam = np.median(chi2) / stats.chi2.isf(0.5, df=1)
    ax.set_xlabel(r"Expected $-\log_{10}$ p")
    ax.set_ylabel(r"Observed $-\log_{10}$ p")
    ax.set_title(f"{title}  (λ={lam:.2f})", fontsize=9)
    ax.set_xlim(0, max(lo, hi) * 1.05)
    ax.set_ylim(0, max(lo, hi) * 1.05)


def savefig(fig, base):
    """Save the figure as both PNG and PDF."""
    fig.savefig(os.path.join(OUT, base + ".png"), dpi=150)
    fig.savefig(os.path.join(OUT, base + ".pdf"))


def main():
    os.makedirs(OUT, exist_ok=True)
    for t in KEY_TRAITS:
        for s in SETS:
            f = os.path.join(SRC, f"{s}_{t}_assoc.csv.gz")
            if not os.path.exists(f):
                continue
            df = pd.read_csv(f)
            df = df[df["p_wald"].notna()]
            df["chr"] = df["chr"].astype(int)
            fig, ax = plt.subplots(1, 2, figsize=(16, 4))
            manhattan(ax[0], df, f"{t} ({s})")
            qq(ax[1], df, f"{t} ({s})")
            fig.tight_layout()
            savefig(fig, f"manhattan_{s}_{t}")
            plt.close(fig)
            print(f"saved {s}_{t}")

    # multi-trait Fisher manhattans
    for s in SETS:
        for block in ["color", "copper"]:
            f = os.path.join(SRC, f"multi_{s}_{block}_fisher_top500.csv")
            if not os.path.exists(f):
                continue
            df = pd.read_csv(f)
            df = df[df["fisher_p"].notna()]
            df["chr"] = df["chr"].astype(int)
            fig, ax = plt.subplots(1, 2, figsize=(16, 4))
            manhattan(ax[0], df, f"{block} multi-Fisher ({s}) top500", pcol="fisher_p")
            p = df["fisher_p"].values
            n = len(p)
            exp = -np.log10(np.linspace(1 / (n + 1), n / (n + 1), n))
            obs = -np.log10(np.sort(p))
            ax[1].scatter(exp, obs, s=4, alpha=0.5, c="#EE6677", linewidths=0)
            lo, hi = max(exp), max(obs)
            ax[1].plot([0, max(lo, hi)], [0, max(lo, hi)], "k--", lw=1)
            ax[1].set_xlabel("Expected")
            ax[1].set_ylabel("Observed")
            ax[1].set_title(f"{block} multi-Fisher ({s}) top500", fontsize=9)
            fig.tight_layout()
            savefig(fig, f"manhattan_multi_{s}_{block}")
            plt.close(fig)
            print(f"saved multi {s} {block}")

    print("ALL DONE ->", OUT)


if __name__ == "__main__":
    main()
