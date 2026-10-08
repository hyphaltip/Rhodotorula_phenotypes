#!/usr/bin/env python3
"""Figures for s1 (species/population), s3 (regimes), s4 (size-matched), s6 (batch) from the report tables."""
import sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, "analysis/carotenoid_stress_vs_baseline/scripts/strat")
from strat_common import *
T = REP / "tables"; F = REP / "figures"; M4 = ["Chromium", "Copper", "Iron", "Lead"]
short = lambda s: s.replace("Rhodotorula ", "R. ")
star = lambda p: "***" if p < .001 else "**" if p < .01 else "*" if p < .05 else ""
# ---- s1 species baseline forest
b = pd.read_csv(T / "s1_species_baseline.csv"); b = b[b.size_adjusted == True]
fig, ax = plt.subplots(1, 4, figsize=(19, 5.2), sharey=False)
for k, m in enumerate(M4):
    s = b[b.Metal == m].sort_values("estimate"); y = np.arange(len(s))
    ax[k].errorbar(s.estimate, y, xerr=1.96 * s.se, fmt="o", color=COL[m], capsize=2); ax[k].axvline(0, color="k", lw=.8)
    ax[k].set_yticks(y); ax[k].set_yticklabels([f"{short(a)} (n={int(n)}) {star(p)}" for a, n, p in zip(s.species, s.n_strains, s.p_BH)], fontsize=8)
    ax[k].set_title(m, fontsize=10); ax[k].set_xlabel("baseline a* difference vs R. mucilaginosa (size-adjusted, 95% CI)")
fig.suptitle("Species effects on baseline a* (mixed model; strain, run, plate random). Stars: BH-adjusted p < .05, .01, .001", fontsize=10)
fig.tight_layout(); fig.savefig(F / "s1_species_baseline_forest.png", dpi=170); plt.close(fig)
# ---- s1 species dose slope forest
d = pd.read_csv(T / "s1_species_dose_slope.csv")
fig, ax = plt.subplots(1, 4, figsize=(19, 5.2))
for k, m in enumerate(M4):
    s = d[d.Metal == m].sort_values("slope_full_range"); y = np.arange(len(s))
    ax[k].errorbar(s.slope_full_range, y, xerr=1.96 * s.se, fmt="o", color=COL[m], capsize=2); ax[k].axvline(0, color="k", lw=.8)
    ax[k].set_yticks(y); ax[k].set_yticklabels([f"{short(a)} (n={int(n)}){' [ref]' if pd.isna(p) else ' ' + star(p)}" for a, n, p in zip(s.species, s.n_strains, s.vs_reference_p)], fontsize=8)
    ax[k].set_title(m, fontsize=10); ax[k].set_xlabel("change in a* over the full dose range, at fixed size (95% CI)")
fig.suptitle("Species-specific dose response of a*. Stars: interaction with R. mucilaginosa, unadjusted p", fontsize=10)
fig.tight_layout(); fig.savefig(F / "s1_species_slope_forest.png", dpi=170); plt.close(fig)
# ---- s1 population
pb = pd.read_csv(T / "s1_population_baseline.csv"); ps = pd.read_csv(T / "s1_population_dose_slope.csv")
fig, ax = plt.subplots(2, 4, figsize=(19, 8))
for k, m in enumerate(M4):
    s = pb[pb.Metal == m]; ax[0, k].errorbar(range(len(s)), s.estimate, yerr=1.96 * s.se, fmt="o", color=COL[m], capsize=2); ax[0, k].axhline(0, color="k", lw=.8)
    ax[0, k].set_xticks(range(len(s))); ax[0, k].set_xticklabels([f"{p.replace('poppop','pop')}\n(n={int(n)})" for p, n in zip(s.population, s.n_strains)], fontsize=8); ax[0, k].set_title(f"{m}: baseline a* vs pop1 (size-adjusted)", fontsize=9)
    s = ps[ps.Metal == m]; ax[1, k].errorbar(range(len(s)), s.slope_full_range, yerr=1.96 * s.se, fmt="o", color=COL[m], capsize=2); ax[1, k].axhline(0, color="k", lw=.8)
    ax[1, k].set_xticks(range(len(s))); ax[1, k].set_xticklabels([f"{p}\n(n={int(n)})" for p, n in zip(s.population, s.n_strains)], fontsize=8); ax[1, k].set_title(f"{m}: a* change over full dose range, by population", fontsize=9)
fig.suptitle("Population structure within R. mucilaginosa (201 strains with a population label)", fontsize=10)
fig.tight_layout(); fig.savefig(F / "s1_population.png", dpi=170); plt.close(fig)
# ---- s3 regimes
r = pd.read_csv(T / "s3_size_ratio_by_dose.csv"); e = pd.read_csv(T / "s3_dose_effects.csv"); sl = pd.read_csv(T / "s3_regime_slopes.csv")
fig, ax = plt.subplots(2, 4, figsize=(19, 8))
for k, m in enumerate(M4):
    x = r[r.Metal == m]; ax[0, k].fill_between(x.conc, x.q25, x.q75, color=COL[m], alpha=.25); ax[0, k].plot(x.conc, x.median_ratio, "-o", color=COL[m]); ax[0, k].axhline(.5, color="k", ls="--", lw=.8)
    for _, z in x.iterrows(): ax[0, k].axvspan(z.conc - .02 * x.conc.max(), z.conc + .02 * x.conc.max(), color=("#999" if z.regime == "inhibitory" else "#fff"), alpha=.25, lw=0)
    ax[0, k].set_title(f"{m}: colony area relative to 0 dose (median, IQR across strains)\ndashed = 0.5 cut; grey = inhibitory", fontsize=8); ax[0, k].set_xlabel("conc"); ax[0, k].set_ylabel("area ratio")
    for adj, ls, lab in ((False, "-", "without size"), (True, "--", "at fixed size")):
        z = e[(e.Metal == m) & (e.size_adjusted == adj)]; ax[1, k].errorbar(z.conc, z.effect_vs_0, yerr=1.96 * z.se, fmt="o" + ls, color=COL[m], alpha=1 if adj else .5, capsize=2, label=lab)
    ax[1, k].axhline(0, color="k", lw=.7); ax[1, k].set_title(f"{m}: a* difference from 0 dose (dose-factor model, 95% CI)", fontsize=8); ax[1, k].set_xlabel("conc"); ax[1, k].legend(fontsize=7)
fig.suptitle("Dose regimes from colony size (threshold 0.5), and a* response by dose", fontsize=10)
fig.tight_layout(); fig.savefig(F / "s3_regimes.png", dpi=170); plt.close(fig)
fig, ax = plt.subplots(figsize=(9, 4.4)); lab = []; i = 0
for m in M4:
    for reg in ("sub-inhibitory", "inhibitory"):
        for j, adj in enumerate((False, True)):
            z = sl[(sl.Metal == m) & (sl.regime == reg) & (sl.size_adjusted == adj)]
            if len(z): ax.bar(i + j * .38, z.slope_per_0_1_of_max_dose.iloc[0] if "slope_per_0_1_of_max_dose" in z else z.iloc[0, 5], .36, yerr=1.96 * z.se.iloc[0], color=("#0072B2", "#D55E00")[j], label=("without size", "at fixed size")[j] if i == 0 else None)
        lab.append(f"{m[:2]}\n{reg[:3]}"); i += 1
ax.set_xticks(np.arange(len(lab)) + .19); ax.set_xticklabels(lab, fontsize=8); ax.axhline(0, color="k", lw=.6); ax.legend(); ax.set_ylabel("a* change per 10% of the metal's max dose"); ax.set_title("Dose slope by regime", fontsize=9)
fig.tight_layout(); fig.savefig(F / "s3_regime_slopes.png", dpi=170); plt.close(fig)
# ---- s4 size-matched
wl = load_wells(species=False); s4 = pd.read_csv(T / "s4_size_matched.csv")
fig, ax = plt.subplots(3, 4, figsize=(19, 12))
for k, m in enumerate(M4):
    g = wl[wl.Metal == m].copy(); q = np.quantile(g.lnA, np.linspace(0, 1, 6)); q[0] -= 1e-6; g["sb"] = pd.cut(g.lnA, q, labels=[f"S{i}" for i in range(1, 6)])
    cov = g.groupby(["sb", "conc"], observed=False).size().unstack(fill_value=0)
    im = ax[0, k].imshow(cov.values, aspect="auto", cmap="Blues"); ax[0, k].set_xticks(range(cov.shape[1])); ax[0, k].set_xticklabels([f"{c:g}" for c in cov.columns], fontsize=7); ax[0, k].set_yticks(range(5)); ax[0, k].set_yticklabels(cov.index)
    for i in range(cov.shape[0]):
        for j in range(cov.shape[1]): ax[0, k].text(j, i, int(cov.values[i, j]), ha="center", va="center", fontsize=6)
    ax[0, k].set_title(f"{m}: wells per size bin x dose (S1 smallest)", fontsize=8)
    for i, sb in enumerate(cov.index):
        z = g[g.sb == sb].groupby("conc").a.agg(["mean", "sem", "size"]); z = z[z["size"] >= 15]
        ax[1, k].errorbar(z.index, z["mean"], yerr=1.96 * z["sem"], fmt="-o", ms=3, color=plt.cm.viridis(i / 4), label=sb, capsize=1)
    ax[1, k].legend(fontsize=7, title="size bin"); ax[1, k].set_xlabel("conc"); ax[1, k].set_ylabel("mean a* (95% CI; bins with >= 15 wells)"); ax[1, k].set_title(f"{m}: a* vs dose within size bins", fontsize=8)
    ax[2, k].boxplot([g[g.sb == f"S{i}"].a for i in range(1, 6)], tick_labels=[f"S{i}" for i in range(1, 6)], showfliers=False, patch_artist=True, boxprops=dict(facecolor=COL[m], alpha=.5))
    ax[2, k].set_title(f"{m}: spread of a* by size bin (small-colony check)", fontsize=8); ax[2, k].set_xlabel("size bin"); ax[2, k].set_ylabel("a*")
fig.suptitle("Size-matched view. Bins are quantiles of ln area within each metal; high doses populate small bins, so coverage is uneven", fontsize=10)
fig.tight_layout(); fig.savefig(F / "s4_size_matched.png", dpi=170); plt.close(fig)
# ---- s6 batch
p = pd.read_csv(T / "s6_per_run_dose_effect.csv"); h = pd.read_csv(T / "s6_run_heterogeneity.csv"); v = pd.read_csv(T / "s6_variance_components.csv")
fig, ax = plt.subplots(1, 5, figsize=(22, 4.4), gridspec_kw=dict(width_ratios=[1, 1, 1, 1, 1.2]))
for k, m in enumerate(M4):
    z = p[p.Metal == m]; ax[k].errorbar(z.dose_effect, range(len(z)), xerr=1.96 * z.se, fmt="o", color=COL[m], capsize=2)
    hh = h[h.Metal == m].iloc[0]; ax[k].axvline(hh.pooled_effect, color="k", ls="--", lw=.8)
    ax[k].set_yticks(range(len(z))); ax[k].set_yticklabels([f"{r} ({int(n)} strains)" for r, n in zip(z.run, z.n_strains)], fontsize=7)
    ax[k].set_title(f"{m}: I2={hh.I2:.2f}, Q p={hh.Q_p:.2g}", fontsize=9); ax[k].set_xlabel("a* change over full dose range, fixed size (95% CI)")
vv = v.set_index("Metal")[["share_strain", "share_run", "share_plate", "share_resid"]].reindex(M4) * 100
vv.plot.bar(stacked=True, ax=ax[4], color=["#0072B2", "#D55E00", "#999999", "#E69F00"]); ax[4].set_ylabel("% of random variance"); ax[4].set_title("Variance components (run = strain subset)", fontsize=9); ax[4].legend(fontsize=7)
fig.suptitle("Batch: dose effect by run (each run is a different strain subset) and variance components", fontsize=10)
fig.tight_layout(); fig.savefig(F / "s6_batch.png", dpi=170); plt.close(fig)
print("figures written")
