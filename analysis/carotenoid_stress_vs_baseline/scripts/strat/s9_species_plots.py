#!/usr/bin/env python3
"""Species-stratified figures: (A) R. mucilaginosa alone (with population), (B) each other species alone against R. mucilaginosa as a grey reference.
Strain is the unit for lines (strain means across replicate wells; error bars = 95% CI across strains). Zinc is not shown (only mucilaginosa has data, one plate per dose)."""
import sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, "analysis/carotenoid_stress_vs_baseline/scripts/strat")
from strat_common import *
F = REP / "figures"; T = REP / "tables"; M4 = ["Chromium", "Copper", "Lead"]; REF = "Rhodotorula mucilaginosa"; MINS = 5
short = lambda s: s.replace("Rhodotorula ", "R. ")
w = load_wells(); w = w[w.Metal.isin(M4) & w.species.notna() & (w.species != "Species Not Found")].copy()
cnt = w.groupby(["Metal", "species"]).strain_id.nunique().unstack(0).fillna(0)
SP = [s for s in cnt.index if s != REF and (cnt.loc[s] >= MINS).sum() >= 2]   # species with >= 5 strains in at least 2 metals (R. graminis: Cr only, 2 doses -> tables only); SP = sorted(SP, key=lambda s: -cnt.loc[s].sum())
PAL = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#000000"]; SC = {s: PAL[i % 8] for i, s in enumerate(SP)}
print("species strata:", {short(s): cnt.loc[s].astype(int).to_dict() for s in SP}); print("mucilaginosa strains:", cnt.loc[REF].astype(int).to_dict())
def strain_means(g):
    return g.groupby(["strain_id", "conc"]).agg(a=("a", "mean"), lnA=("lnA", "mean")).reset_index()
def line(g, col):
    s = strain_means(g); z = s.groupby("conc")[col].agg(["mean", "sem", "size"]); return z[z["size"] >= 3]
# ---------- A1: R. mucilaginosa a* vs size by dose
mu = w[w.species == REF]
fig, ax = plt.subplots(1, 3, figsize=(19, 4.4), sharey=True)
for k, m in enumerate(M4):
    g = mu[mu.Metal == m]; doses = sorted(g.conc.unique()); edges = np.quantile(g.lnA, np.linspace(0, 1, 11)); edges[-1] += 1e-9
    for i, d in enumerate(doses):
        x = g[g.conc == d]; b = np.digitize(x.lnA, edges) - 1; s = x.groupby(b).agg(l=("lnA", "median"), a=("a", "median"), n=("a", "size")); s = s[s.n >= 15]
        ax[k].plot(s.l, s.a, "-o", ms=3, color=plt.cm.viridis(i / max(len(doses) - 1, 1)), label=f"{d:g}")
    ax[k].set_title(f"{m}: R. mucilaginosa, {g.strain_id.nunique()} strains, {len(g):,} wells", fontsize=9); ax[k].set_xlabel("ln colony area (px), well median"); ax[k].legend(title="conc", fontsize=7, title_fontsize=7)
ax[0].set_ylabel("a* (well median)"); fig.suptitle("R. mucilaginosa only: a* against colony size by dose", fontsize=10); fig.tight_layout(); fig.savefig(F / "sp_muc_a_vs_size.png", dpi=170); plt.close(fig)
# ---------- B1/B2: every species, dose response of a* and ln area, mucilaginosa grey reference
for col, fn, yl in (("a", "sp_species_dose_response_a", "mean a* (strain means, 95% CI)"), ("lnA", "sp_species_dose_response_size", "mean ln colony area")):
    fig, ax = plt.subplots(len(SP), 3, figsize=(17, 2.6 * len(SP)), sharex="col")
    for k, m in enumerate(M4):
        ref = line(w[(w.Metal == m) & (w.species == REF)], col)
        for r_, s in enumerate(SP):
            a_ = ax[r_, k]; a_.plot(ref.index, ref["mean"], "-", color="#999999", lw=1.4, label="R. mucilaginosa")
            g = w[(w.Metal == m) & (w.species == s)]; n = g.strain_id.nunique()
            if n >= MINS:
                z = line(g, col); a_.errorbar(z.index, z["mean"], yerr=1.96 * z["sem"], fmt="-o", ms=3, color=SC[s], capsize=1.5)
            else: a_.text(.5, .5, f"n = {n} strains (< {MINS})", transform=a_.transAxes, ha="center", fontsize=8, color="#777")
            a_.set_title(f"{short(s)}, {m} (n={n})", fontsize=8)
            if k == 0: a_.set_ylabel(yl, fontsize=7)
        ax[-1, k].set_xlabel("conc")
    fig.suptitle(("a*" if col == "a" else "Colony size (ln area)") + " by dose, each species (colour) against R. mucilaginosa (grey)", fontsize=10); fig.tight_layout(rect=[0, 0, 1, 0.975]); fig.savefig(F / (fn + ".png"), dpi=150); plt.close(fig)
# ---------- B3: model-based effect at fixed size by stratum
dfc = pd.read_csv(T / "s9_stratum_dose_factor.csv")
fig, ax = plt.subplots(len(SP), 3, figsize=(17, 2.6 * len(SP)), sharex="col")
for k, m in enumerate(M4):
    ref = dfc[(dfc.Metal == m) & (dfc.stratum == REF)]
    for r_, s in enumerate(SP):
        a_ = ax[r_, k]; a_.plot(ref.conc, ref.effect_at_fixed_size, "-", color="#999999", lw=1.4); a_.axhline(0, color="k", lw=.5)
        z = dfc[(dfc.Metal == m) & (dfc.stratum == s)]
        if len(z): a_.errorbar(z.conc, z.effect_at_fixed_size, yerr=1.96 * z.se, fmt="-o", ms=3, color=SC[s], capsize=1.5)
        else: a_.text(.5, .5, "not fitted (< 5 strains)", transform=a_.transAxes, ha="center", fontsize=8, color="#777")
        a_.set_title(f"{short(s)}, {m}", fontsize=8)
        if k == 0: a_.set_ylabel("a* vs 0 dose, fixed size", fontsize=7)
    ax[-1, k].set_xlabel("conc")
fig.suptitle("Mixed-model dose effect on a* at fixed size, each species (colour) against R. mucilaginosa (grey); 95% CI", fontsize=10); fig.tight_layout(rect=[0, 0, 1, 0.975]); fig.savefig(F / "sp_species_effect_fixed_size.png", dpi=150); plt.close(fig)
# ---------- B4: a* vs ln area scatter by species, one figure per metal
ALL = [REF] + SP
for m in M4:
    nr = -(-len(ALL) // 3); fig, ax = plt.subplots(nr, 3, figsize=(15, 4 * nr), sharex=True, sharey=True, squeeze=False); g = w[w.Metal == m]
    for a_ in ax.flat[len(ALL):]: a_.axis('off')
    for i, s in enumerate(ALL):
        a_ = ax[i // 3, i % 3]; x = g[g.species == s]; n = x.strain_id.nunique()
        if n >= MINS:
            sc = a_.scatter(x.lnA, x.a, c=x.dose_s, cmap="viridis", s=7, alpha=.6, vmin=0, vmax=1)
        a_.set_title(f"{short(s)}: {n} strains, {len(x):,} wells" if n >= MINS else f"{short(s)}: {n} strains (< {MINS})", fontsize=9)
    for a_ in ax[nr - 1]: a_.set_xlabel("ln colony area (px)")
    for a_ in ax[:, 0]: a_.set_ylabel("a* (well median)")
    fig.colorbar(sc, ax=ax, label="dose (fraction of max)", shrink=.6); fig.suptitle(f"{m}: a* against colony size, by species (colour = dose)", fontsize=11); fig.savefig(F / f"sp_{m.lower()}_a_vs_size_by_species.png", dpi=140); plt.close(fig)
# ---------- B5: baseline vs top dose per strain, by species
fig, ax = plt.subplots(len(SP), 3, figsize=(17, 2.9 * len(SP)))
for k, m in enumerate(M4):
    g0 = w[w.Metal == m]; dm = g0.conc.max()
    for r_, s in enumerate(SP):
        a_ = ax[r_, k]; ref = strain_means(g0[g0.species == REF]); ref = ref[ref.conc.isin([0, dm])].pivot(index="strain_id", columns="conc", values="a").dropna()
        a_.scatter(ref[0.0], ref[dm], s=5, color="#cccccc"); gs = g0[g0.species == s]
        t = strain_means(gs); t = t[t.conc.isin([0, dm])].pivot(index="strain_id", columns="conc", values="a").dropna() if len(gs) else pd.DataFrame()
        if len(t) >= 3: a_.scatter(t[0.0], t[dm], s=18, color=SC[s])
        lim = [ref.min().min(), ref.max().max()]; a_.plot(lim, lim, "k--", lw=.7); a_.set_title(f"{short(s)}, {m} (n={len(t)})", fontsize=8)
        if k == 0: a_.set_ylabel(f"a*, top dose", fontsize=7)
    ax[-1, k].set_xlabel("a*, 0 dose")
fig.suptitle("Baseline against top-dose a* per strain, each species (colour) over R. mucilaginosa (grey); dashed = no change", fontsize=10); fig.tight_layout(rect=[0, 0, 1, 0.975]); fig.savefig(F / "sp_species_baseline_vs_top.png", dpi=150); plt.close(fig)
# ---------- table: strain counts
cnt.loc[[REF] + SP].astype(int).reset_index().rename(columns={"species": "species"}).to_csv(T / "s9_species_strain_counts.csv", index=False)
print("done")
