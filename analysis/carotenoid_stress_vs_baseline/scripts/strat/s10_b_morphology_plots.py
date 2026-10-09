#!/usr/bin/env python3
"""b* and morphology figures from wells_traits.csv; strain is the unit for lines (strain means, 95% CI across strains)."""
import sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy import stats
sys.path.insert(0, "analysis/carotenoid_stress_vs_baseline/scripts/strat")
from strat_common import *
F = REP / "figures"; T = REP / "tables"; M4 = ["Chromium", "Copper", "Lead"]; REF = "Rhodotorula mucilaginosa"
LAB = {"L": "L* (lightness)", "a": "a* (red-green)", "b": "b* (yellow-blue)", "chroma": "chroma sqrt(a*^2+b*^2)", "hue_deg": "hue angle (deg; 0 = red, 90 = yellow)", "sat": "HSV saturation", "val": "HSV value",
       "circ": "circularity", "solid": "solidity", "ecc": "eccentricity", "compact": "compactness", "extent": "extent", "aspect": "aspect ratio (major/minor)"}
COLT = ["L", "a", "b", "chroma", "hue_deg"]; MORPH = ["circ", "solid", "ecc", "compact", "extent", "aspect"]
w = pd.read_csv(R / "wells_traits.csv"); st = pd.read_csv(R / "strat/strain_table.csv", dtype={"strain_id": str}); w["strain_id"] = w.strain_id.astype(str)
w = w.merge(st[["strain_id", "species"]], on="strain_id", how="left"); w = w[w.Metal.isin(M4)].copy()
w["grp"] = np.where(w.species == REF, "R. mucilaginosa", np.where(w.species.notna() & (w.species != "Species Not Found"), "other species", "no species"))
print("wells:", len(w), w.grp.value_counts().to_dict())
def line(g, col):
    s = g.groupby(["strain_id", "conc"])[col].mean().reset_index(); z = s.groupby("conc")[col].agg(["mean", "sem", "size"]); return z[z["size"] >= 3]
GC = {"R. mucilaginosa": "#0072B2", "other species": "#D55E00"}
def grid(traits, fn, title):
    fig, ax = plt.subplots(len(traits), 3, figsize=(17, 2.7 * len(traits)), sharex="col")
    for k, m in enumerate(M4):
        g = w[w.Metal == m]
        for r_, tr in enumerate(traits):
            for gn in ("R. mucilaginosa", "other species"):
                z = line(g[g.grp == gn], tr)
                if len(z): ax[r_, k].errorbar(z.index, z["mean"], yerr=1.96 * z["sem"], fmt="-o", ms=3, color=GC[gn], capsize=1.5, label=f"{gn} ({g[g.grp == gn].strain_id.nunique()} strains)")
            if k == 0: ax[r_, k].set_ylabel(LAB[tr], fontsize=7)
            ax[r_, k].set_title(f"{m}: {tr}", fontsize=8)
            if r_ == 0: ax[r_, k].legend(fontsize=6)
        ax[-1, k].set_xlabel("conc")
    fig.suptitle(title, fontsize=10); fig.tight_layout(); fig.savefig(F / fn, dpi=150); plt.close(fig)
grid(COLT, "bm_colour_response.png", "Colour traits by dose (strain means, 95% CI across strains): R. mucilaginosa and all other species")
grid(MORPH, "bm_morphology_response.png", "Morphology by dose (strain means, 95% CI across strains): R. mucilaginosa and all other species. Shape metrics of very small colonies are noisy")
# a*-b* plane
fig, ax = plt.subplots(1, 3, figsize=(19, 4.8))
for k, m in enumerate(M4):
    g = w[w.Metal == m]; sm = g.groupby(["strain_id", "conc"])[["a", "b"]].mean().reset_index(); doses = sorted(g.conc.unique())
    for i, d in enumerate(doses):
        z = sm[sm.conc == d]; ax[k].scatter(z.a, z.b, s=6, alpha=.25, color=plt.cm.viridis(i / max(len(doses) - 1, 1)))
    mm = sm.groupby("conc")[["a", "b"]].mean(); ax[k].plot(mm.a, mm.b, "k-", lw=1.2)
    for i, (d, r) in enumerate(mm.iterrows()): ax[k].scatter(r.a, r.b, s=55, color=plt.cm.viridis(i / max(len(doses) - 1, 1)), edgecolor="k", zorder=3); ax[k].annotate(f"{d:g}", (r.a, r.b), fontsize=7, xytext=(4, 4), textcoords="offset points")
    ax[k].set_xlabel("a* (red-green)"); ax[k].set_ylabel("b* (yellow-blue)"); ax[k].set_title(f"{m}: strain means by dose; black line = mean over strains", fontsize=8)
fig.suptitle("Colour-plane trajectory under stress (colour = dose, labels = conc)", fontsize=10); fig.tight_layout(); fig.savefig(F / "bm_ab_plane.png", dpi=170); plt.close(fig)
# strain-level correlation heatmaps: morphology vs colour at dose 0, and change at top dose
MC = MORPH + ["lnA"]; CC = ["a", "b", "L", "chroma", "hue_deg"]
def cormat(df, rowt, colt):
    return pd.DataFrame({c: {r: stats.spearmanr(df[r], df[c], nan_policy="omit")[0] for r in rowt} for c in colt})
fig, ax = plt.subplots(2, 3, figsize=(18, 11.5), gridspec_kw=dict(hspace=0.55, wspace=0.4)); rows = []   # room for rotated x labels between the rows
for k, m in enumerate(M4):
    g = w[w.Metal == m]; dm = g.conc.max(); s0 = g[g.conc == 0].groupby("strain_id")[MC + CC].mean(); st_ = g[g.conc == dm].groupby("strain_id")[MC + CC].mean()
    dl = (st_ - s0).dropna()
    for r_, (df, ttl) in enumerate(((s0.dropna(), "baseline (0 dose), strain means"), (dl, f"change (dose {dm:g} minus 0)"))):
        cm = cormat(df, MC, CC); im = ax[r_, k].imshow(cm.values, cmap="RdBu_r", vmin=-1, vmax=1); ax[r_, k].set_xticks(range(len(CC))); ax[r_, k].set_xticklabels(CC, rotation=45, ha="right", fontsize=8)
        ax[r_, k].set_yticks(range(len(MC))); ax[r_, k].set_yticklabels(MC, fontsize=8); ax[r_, k].set_title(f"{m}\n{ttl}; n={len(df)}", fontsize=8)
        for i in range(cm.shape[0]):
            for j in range(cm.shape[1]): ax[r_, k].text(j, i, f"{cm.values[i, j]:.2f}", ha="center", va="center", fontsize=6)
        cm["metal"] = m; cm["scope"] = "baseline" if r_ == 0 else "change"; cm["morph"] = cm.index; rows.append(cm)
fig.suptitle("Strain-level correlation of morphology (rows) with colour (columns)", fontsize=10); fig.subplots_adjust(top=0.9, bottom=0.08, left=0.06, right=0.9); fig.colorbar(im, cax=fig.add_axes([0.925, 0.3, 0.012, 0.4]), label="Spearman rho"); fig.savefig(F / "bm_morph_colour_corr.png", dpi=150); plt.close(fig)
pd.concat(rows).to_csv(T / "s10_morph_colour_corr.csv", index=False)
# baseline b* and circularity by species
sp = w[(w.grp != "no species") & w.species.notna()]; cn = sp.groupby(["Metal", "species"]).strain_id.nunique().unstack(0).fillna(0)
fig, ax = plt.subplots(2, 3, figsize=(19, 9))
for k, m in enumerate(M4):
    g = sp[(sp.Metal == m) & (sp.conc == 0)]; sm = g.groupby("strain_id").agg(b=("b", "mean"), circ=("circ", "mean"), species=("species", "first"))
    keep = [s for s in cn.index if cn.loc[s, m] >= 5]
    for r_, col in enumerate(("b", "circ")):
        order = sm[sm.species.isin(keep)].groupby("species")[col].median().sort_values().index
        ax[r_, k].boxplot([sm[sm.species == o][col] for o in order], tick_labels=[f"{o.replace('Rhodotorula ', 'R. ')} ({(sm.species == o).sum()})" for o in order], vert=False, showfliers=False)
        ax[r_, k].set_title(f"{m}: baseline {LAB[col]}", fontsize=8); ax[r_, k].tick_params(axis="y", labelsize=7)
fig.suptitle("Baseline b* and circularity by species (strain means at 0 dose; species with >= 5 strains)", fontsize=10); fig.tight_layout(); fig.savefig(F / "bm_species_baseline_b_circ.png", dpi=150); plt.close(fig)
# model-based dose effects heatmap (needs s10_trait_dose_effects.csv)
try:
    e = pd.read_csv(T / "s10_trait_dose_effects.csv")
    fig, ax = plt.subplots(2, 2, figsize=(17, 14), gridspec_kw=dict(hspace=0.4, wspace=0.18)); star = lambda p: "***" if p < .001 else "**" if p < .01 else "*" if p < .05 else ""
    for r_, ds in enumerate(("main: largest object >= 2000 px", "no area filter")):
        for c_, mod in enumerate(("total", "at fixed size")):
            z = e[(e.dataset == ds) & (e.model == mod)]; pv = z.pivot(index="trait", columns="Metal", values="effect_sd").reindex(COLT + ["sat", "val"] + MORPH)[M4]; pp = z.pivot(index="trait", columns="Metal", values="p").reindex(pv.index)[M4]
            im = ax[r_, c_].imshow(pv.values, cmap="RdBu_r", vmin=-3, vmax=3, aspect="auto"); ax[r_, c_].set_xticks(range(3)); ax[r_, c_].set_xticklabels(M4, rotation=45, ha="right"); ax[r_, c_].set_yticks(range(len(pv))); ax[r_, c_].set_yticklabels(pv.index)
            for i in range(pv.shape[0]):
                for j in range(3):
                    v = pv.values[i, j]
                    if not np.isnan(v): ax[r_, c_].text(j, i, f"{v:.1f}{star(pp.values[i, j])}", ha="center", va="center", fontsize=9)
            ax[r_, c_].set_title(f"{ds}; {mod}", fontsize=9)
    fig.subplots_adjust(top=0.93, bottom=0.08, left=0.07, right=0.9); fig.colorbar(im, cax=fig.add_axes([0.92, 0.3, 0.014, 0.4]), label="dose effect (SD of 0-dose wells)"); fig.suptitle("Mixed-model dose effect (0 to max dose) on colour and morphology traits, in SD units. Stars: p < .05, .01, .001 (unadjusted)", fontsize=10)
    fig.savefig(F / "bm_dose_effect_heatmap.png", dpi=150); plt.close(fig)
except FileNotFoundError: print("no dose-effect table yet")
print("done")
