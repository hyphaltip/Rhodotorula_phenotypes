#!/usr/bin/env python3
"""Dose-0 consistency: figures and tables (strain-mean agreement between screens, plate/run batch, within-strain heterogeneity)."""
import itertools
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy import stats
R = Path("analysis/gwas_dh4148/results"); T = Path("analysis/gwas_dh4148/report/tables"); F = Path("analysis/gwas_dh4148/report/figures")
w = pd.read_csv(R / "dose0_wells.csv.gz"); st = pd.read_csv(R / "strain_traits.csv")[["strain_id", "sample_name", "species", "ploidy_status", "gwas_panel"]].drop_duplicates("strain_id")
w = w.merge(st, on="strain_id", how="left")
M3 = ["Chromium", "Copper", "Lead"]; COL = {"Chromium": "#0072B2", "Copper": "#D55E00", "Lead": "#CC79A7", "Iron": "#009E73", "Zinc": "#E69F00"}
TR = {"a": "a* (late window)", "lnA": "ln colony area (late window)", "rgr": "max growth rate (1/h)"}
log = lambda m: print(m, flush=True)
# 1. variance components
vc = pd.read_csv(T / "dose0_variance_components.csv"); v = vc[vc.scope == "Cr+Cu+Pb pooled"].set_index("trait")
fig, ax = plt.subplots(1, 2, figsize=(12, 4.2))
parts = ["share_strain", "share_metal", "share_run", "share_plate", "share_resid"]; cols = ["#0072B2", "#E69F00", "#009E73", "#999999", "#D55E00"]
(v.loc[list(TR), parts] * 100).rename(columns=lambda c: c.replace("share_", "")).plot.bar(stacked=True, ax=ax[0], color=cols); ax[0].set_ylabel("% of dose-0 variance"); ax[0].set_title("Pooled Cr + Cu + Pb dose-0 wells", fontsize=10); ax[0].set_xticklabels([TR[t] for t in TR], rotation=15, ha="right", fontsize=8)
s = vc[vc.scope.isin(M3 + ["Cr+Cu+Pb pooled"]) & (vc.trait == "a")].set_index("scope").loc[M3 + ["Cr+Cu+Pb pooled"]]
(s[parts[:1] + parts[2:]] * 100).rename(columns=lambda c: c.replace("share_", "")).plot.bar(stacked=True, ax=ax[1], color=[cols[0]] + cols[2:]); ax[1].set_title("a* by screen", fontsize=10); ax[1].set_ylabel("% of dose-0 variance"); ax[1].tick_params(axis="x", rotation=15, labelsize=8)
fig.tight_layout(); fig.savefig(F / "d0_variance_components.png", dpi=170); plt.close(fig)
# 2. strain means across screens
sm = {}
for tr in TR: sm[tr] = w[w.Metal.isin(M3)].groupby(["strain_id", "Metal"])[tr].mean().unstack()
rows = []
fig, ax = plt.subplots(3, 3, figsize=(12, 11))
for i, tr in enumerate(TR):
    for j, (m1, m2) in enumerate(itertools.combinations(M3, 2)):
        z = sm[tr][[m1, m2]].dropna(); r, p = stats.spearmanr(z[m1], z[m2]); rows.append(dict(trait=tr, screen_1=m1, screen_2=m2, n_strains=len(z), spearman=r, p=p, pearson=np.corrcoef(z[m1], z[m2])[0, 1]))
        ax[i, j].scatter(z[m1], z[m2], s=8, alpha=.5, color="#0072B2"); lim = [min(z.min()), max(z.max())]; ax[i, j].plot(lim, lim, "k--", lw=.7)
        ax[i, j].set_xlabel(f"{m1} dose-0 strain mean"); ax[i, j].set_ylabel(f"{m2}"); ax[i, j].set_title(f"{TR[tr]}: rho = {r:.2f} (n = {len(z)})", fontsize=8)
fig.tight_layout(); fig.savefig(F / "d0_strain_means_between_screens.png", dpi=170); plt.close(fig)
pd.DataFrame(rows).to_csv(T / "dose0_between_screen_correlations.csv", index=False)
# 3. plate / run batch
p = w[w.Metal.isin(M3 + ["Iron", "Zinc"])].groupby(["Metal", "run_number", "plate_position"]).agg(a=("a", "mean"), lnA=("lnA", "mean"), rgr=("rgr", "mean"), n=("a", "size")).reset_index()
p.to_csv(T / "dose0_plate_means.csv", index=False)
fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
for k, tr in enumerate(TR):
    for m, g in p.groupby("Metal"):
        x = [f"{m[:2]}\n{r[-3:]}" for r in g.run_number]; ax[k].scatter(np.arange(len(g)) + list(COL).index(m) * 0.0, g[tr], s=14, color=COL[m], label=m if k == 0 else None)
    ax[k].set_title(TR[tr] + ": mean of each dose-0 plate", fontsize=9); ax[k].set_xlabel("plates, ordered by screen and run")
fig.legend(*ax[0].get_legend_handles_labels(), loc="lower center", ncol=5, fontsize=8); fig.tight_layout(rect=[0, .06, 1, 1]); fig.savefig(F / "d0_plate_means.png", dpi=170); plt.close(fig)
# 4. strain heterogeneity: within-strain SD of dose-0 wells
h = w[w.Metal.isin(M3)].groupby("strain_id").agg(n_wells=("a", "size"), a_mean=("a", "mean"), a_sd=("a", "std"), lnA_mean=("lnA", "mean"), lnA_sd=("lnA", "std"), rgr_mean=("rgr", "mean"), rgr_sd=("rgr", "std")).reset_index().merge(st, on="strain_id")
h = h[h.n_wells >= 6]; h["a_cv"] = h.a_sd / h.a_mean; h["rgr_cv"] = h.rgr_sd / h.rgr_mean
h["a_sd_z"] = (h.a_sd - h.a_sd.median()) / (1.4826 * (h.a_sd - h.a_sd.median()).abs().median())
h.sort_values("a_sd_z", ascending=False).to_csv(T / "dose0_strain_heterogeneity.csv", index=False)
fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
ax[0].scatter(h.a_mean, h.a_sd, s=9, alpha=.5, c=np.where(h.gwas_panel, "#D55E00", "#0072B2")); ax[0].set_xlabel("strain mean a* at dose 0"); ax[0].set_ylabel("SD of a* across its dose-0 wells"); ax[0].set_title("Within-strain spread (red = GWAS panel)", fontsize=9)
for k, (c, lab) in enumerate((("lnA_sd", "ln area"), ("rgr_sd", "growth rate"))): ax[k + 1].scatter(h[c.replace("_sd", "_mean")], h[c], s=9, alpha=.5, c=np.where(h.gwas_panel, "#D55E00", "#0072B2")); ax[k + 1].set_xlabel(f"strain mean {lab}"); ax[k + 1].set_ylabel(f"SD of {lab}")
fig.tight_layout(); fig.savefig(F / "d0_strain_heterogeneity.png", dpi=170); plt.close(fig)
log(h.sort_values("a_sd_z", ascending=False).head(12)[["sample_name", "species", "n_wells", "a_mean", "a_sd", "a_sd_z"]].round(2).to_string(index=False))
log(f"strains with >= 6 dose-0 wells: {len(h)}; SD of a* across wells: median {h.a_sd.median():.2f}; strains with robust z > 3: {(h.a_sd_z > 3).sum()}")
