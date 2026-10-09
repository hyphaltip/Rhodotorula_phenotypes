#!/usr/bin/env python3
"""Figures from the well table and the mixed-model outputs (replicates carried through as error bars)."""
import argparse, sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
METALS = ["Chromium", "Copper", "Lead"]
BLUE, ORANGE = "#0072B2", "#D55E00"

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out", default="analysis/carotenoid_stress_vs_baseline/results"); a = ap.parse_args()
    out = Path(a.out); fig_dir = out / "figures"; fig_dir.mkdir(exist_ok=True)
    w = pd.read_csv(out / "wells.csv"); bl = pd.read_csv(out / "strain_blups.csv"); ms = pd.read_csv(out / "mixed_model_summary.csv")
    # Fig 1: a* vs ln area by dose; each point = median over replicate wells of one strain x dose
    fig, ax = plt.subplots(1, 3, figsize=(21, 4.3), sharey=True)
    for k, m in enumerate(METALS):
        g = w[w.Metal == m]; doses = sorted(g.conc.unique())
        edges = np.quantile(g.lnA, np.linspace(0, 1, 11)); edges[-1] += 1e-9
        for i, d in enumerate(doses):
            x = g[g.conc == d]; b = np.digitize(x.lnA, edges) - 1
            s = x.groupby(b).agg(l=("lnA", "median"), a=("a", "median"), n=("a", "size")); s = s[s.n >= 15]
            ax[k].plot(s.l, s.a, "-o", ms=3, color=plt.cm.viridis(i / max(len(doses) - 1, 1)), label=f"{d:g}")
        ax[k].set_title(f"{m}: {len(g):,} wells", fontsize=10); ax[k].set_xlabel("ln colony area (px), well median")
        ax[k].legend(title="conc (unit not given)", fontsize=7, title_fontsize=7)
    ax[0].set_ylabel("a* (well median)")
    fig.suptitle("a* vs colony size by dose. Overlapping curves = size effect only; separated curves = a* shift at the same size", fontsize=10)
    fig.tight_layout(); fig.savefig(fig_dir / "fig1_astar_vs_size_by_dose.png", dpi=170); plt.close(fig)
    # Fig 2: strain baseline vs top dose, mean of replicate wells +- SE (both axes)
    fig, ax = plt.subplots(2, 3, figsize=(21, 8.3))
    for k, m in enumerate(METALS):
        g = w[w.Metal == m]
        n0 = set(g[g.conc == 0].strain_id)
        dmax = max(d for d in g.conc.unique() if d > 0 and len(n0 & set(g[g.conc == d].strain_id)) >= 30)  # highest dose sharing >= 30 strains with 0 dose
        for r_, (col, lab) in enumerate((("a", "a*"), ("lnA", "ln area"))):
            t0 = g[g.conc.isin([0, dmax])].groupby(["strain_id", "conc"])[col].agg(["mean", "sem", "size"]).unstack("conc")
            t = t0.dropna(subset=[("mean", 0.0), ("mean", dmax)])
            print(f"fig2 {m} {col}: strains with both doses {len(t)} of {len(t0)}; SE undefined (1 well) for {int(t[('sem', 0.0)].isna().sum())} at 0 dose", flush=True)
            if len(t) < 5:
                ax[r_, k].set_title(f"{m}: too few strains with both doses ({len(t)})", fontsize=8); continue
            x, y = t[("mean", 0.0)], t[("mean", dmax)]
            ax[r_, k].errorbar(x, y, xerr=t[("sem", 0.0)].fillna(0), yerr=t[("sem", dmax)].fillna(0), fmt="o", ms=3, alpha=.6, color=BLUE, ecolor="#9ab", lw=.6)
            lim = [min(x.min(), y.min()), max(x.max(), y.max())]; ax[r_, k].plot(lim, lim, "k--", lw=1)
            rho = stats.spearmanr(x, y)[0]
            ax[r_, k].set_xlabel(f"{lab}, 0 dose"); ax[r_, k].set_ylabel(f"{lab}, dose {dmax:g}")
            ax[r_, k].set_title(f"{m}: {len(x)} strains, rho={rho:.2f}, median change {(y-x).median():+.2f}\n(mean of replicate wells, bars = SE; wells/strain/dose median {t[('size', 0.0)].median():.0f})", fontsize=8)
    fig.suptitle("Per-strain a* (top) and colony size (bottom): unstressed vs highest dose. Dashed = no change", fontsize=10)
    fig.tight_layout(); fig.savefig(fig_dir / "fig2_baseline_vs_stressed.png", dpi=170); plt.close(fig)
    # Fig 3: model-based strain effects, total (M0) and at fixed colony size (M1)
    fig, ax = plt.subplots(2, 3, figsize=(21, 8.2))
    for k, m in enumerate(METALS):
        for r_, mod in enumerate(("M0", "M1")):
            b = bl[(bl.Metal == m) & (bl.model == mod)]; r = ms[(ms.Metal == m) & (ms.model == mod)].iloc[0]
            if r.singular:
                ax[r_, k].set_title(f"{m} {mod}: SINGULAR fit, strain-slope BLUPs not interpretable", fontsize=8); continue
            ax[r_, k].scatter(b.blup_intercept, b.blup_dose_slope, s=14, alpha=.7, color=ORANGE)
            ax[r_, k].axhline(0, color="k", lw=.6); ax[r_, k].axvline(0, color="k", lw=.6)
            ax[r_, k].set_xlabel("strain baseline a* (deviation from mean)" + (", at fixed size" if mod == "M1" else ""))
            ax[r_, k].set_ylabel("strain-specific extra change in a* at top dose" + (" (fixed size)" if mod == "M1" else ""), fontsize=8)
            ax[r_, k].set_title(f"{m} {mod}: corr(baseline, induction)={r.cor_intercept_slope:.2f}; n strains={len(b)}", fontsize=8)
    fig.suptitle("Model-based strain effects. M0 = total effect of dose; M1 = effect at the same colony size (size is partly an effect of stress, so M1 is a direct effect, not the whole induction)", fontsize=9)
    fig.tight_layout(); fig.savefig(fig_dir / "fig3_strain_baseline_vs_induction.png", dpi=170); plt.close(fig)
    # Fig 4: share of non-fixed variance due to strain, at 0 dose and at top dose (random slope makes this dose-dependent)
    v = ms[ms.model == "M1"].set_index("Metal")[["repeatability_dose0", "repeatability_dose1"]].reindex(METALS).dropna()
    fig, ax = plt.subplots(figsize=(7.5, 4.2)); v.plot.bar(ax=ax, color=[BLUE, ORANGE]); ax.set_ylabel("strain variance / (strain + plate + residual)")
    ax.legend(["at 0 dose", "at top dose"]); ax.set_title("Repeatability of strain a* across replicate wells (size-adjusted model)", fontsize=9)
    fig.tight_layout(); fig.savefig(fig_dir / "fig4_repeatability.png", dpi=170); plt.close(fig)
    print("figures written")

if __name__ == "__main__":
    sys.exit(main())
