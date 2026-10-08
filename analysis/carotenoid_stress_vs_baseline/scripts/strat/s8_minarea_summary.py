#!/usr/bin/env python3
"""Compare dose effects on a* across minimum-colony-area thresholds (0 = no filter, plus 1000/2000/3000/5000 px)."""
import sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, "analysis/carotenoid_stress_vs_baseline/scripts/strat")
from strat_common import *
TH = [0, 1000, 2000, 3000, 5000]; M4 = ["Chromium", "Copper", "Iron", "Lead"]
mm, m2, wl, topd = [], [], [], []
for A in TH:
    d = R if A == 0 else R / f"minarea_{A}"; wf = R / ("wells.csv" if A == 0 else f"wells_min{A}.csv")
    x = pd.read_csv(d / "mixed_model_summary.csv"); x["min_area"] = A; mm.append(x[x.Metal.isin(M4)])
    y = pd.read_csv(d / "dose_factor_effects_M2.csv"); y["min_area"] = A; m2.append(y[y.Metal.isin(M4)])
    w = pd.read_csv(wf); w = w[w.Metal.isin(M4)]
    c = w.groupby(["Metal", "conc"]).size().reset_index(name="wells"); c["min_area"] = A; wl.append(c)
    for m in M4:
        g = w[w.Metal == m]; dm = g.conc.max(); t = g[g.conc == dm].groupby("strain_id").a.mean()
        topd.append(dict(min_area=A, Metal=m, top_dose=dm, n_strains_top=len(t), top_dose_mean_a=t.mean(), top_dose_sd_across_strains=t.std(),
                         median_area_top_px=np.exp(g[g.conc == dm].lnA).median()))
mm = pd.concat(mm); m2 = pd.concat(m2); wl = pd.concat(wl); topd = pd.DataFrame(topd)
mm.to_csv(REP / "tables/s8_minarea_models.csv", index=False); m2.to_csv(REP / "tables/s8_minarea_dose_factor.csv", index=False); wl.to_csv(REP / "tables/s8_minarea_wells_by_dose.csv", index=False); topd.to_csv(REP / "tables/s8_minarea_topdose.csv", index=False)
fig, ax = plt.subplots(3, 4, figsize=(19, 12))
for k, m in enumerate(M4):
    z = mm[(mm.Metal == m) & (mm.model == "M1")].sort_values("min_area"); z0 = mm[(mm.Metal == m) & (mm.model == "M0")].sort_values("min_area")
    ax[0, k].errorbar(z.min_area, z.dose_effect, yerr=1.96 * z.dose_se, fmt="o-", color=COL[m], label="at fixed size", capsize=2); ax[0, k].errorbar(z0.min_area, z0.dose_effect, yerr=1.96 * z0.dose_se, fmt="s--", color=COL[m], alpha=.5, label="total", capsize=2)
    ax[0, k].axhline(0, color="k", lw=.6); ax[0, k].set_xlabel("minimum colony area (px)"); ax[0, k].set_ylabel("a* change over full dose range"); ax[0, k].set_title(f"{m}: dose effect vs area threshold", fontsize=9); ax[0, k].legend(fontsize=7)
    for A, col in zip(TH, plt.cm.viridis(np.linspace(0, .9, len(TH)))):
        q = m2[(m2.Metal == m) & (m2.min_area == A)]
        ax[1, k].errorbar(q.conc, q.effect_at_fixed_size, yerr=1.96 * q.se, fmt="o-", ms=3, color=col, label=f">= {A} px" if A else "no filter", capsize=1)
    ax[1, k].axhline(0, color="k", lw=.6); ax[1, k].set_xlabel("conc"); ax[1, k].set_ylabel("a* difference from 0 dose (fixed size)"); ax[1, k].legend(fontsize=7); ax[1, k].set_title(f"{m}: per-dose effect by threshold", fontsize=9)
    q = wl[wl.Metal == m].pivot(index="conc", columns="min_area", values="wells").fillna(0)
    for A, col in zip(TH, plt.cm.viridis(np.linspace(0, .9, len(TH)))): ax[2, k].plot(q.index, q[A], "o-", ms=3, color=col, label=f">= {A} px" if A else "no filter")
    ax[2, k].set_yscale("symlog"); ax[2, k].set_xlabel("conc"); ax[2, k].set_ylabel("wells retained"); ax[2, k].legend(fontsize=7); ax[2, k].set_title(f"{m}: wells retained by dose", fontsize=9)
fig.suptitle("Sensitivity to a minimum colony area (well-images whose largest object is smaller are dropped)", fontsize=10)
fig.tight_layout(); fig.savefig(REP / "figures/s8_minarea.png", dpi=170); plt.close(fig)
print(mm[mm.model == "M1"][["Metal", "min_area", "n_wells", "dose_effect", "dose_se", "dose_p", "singular"]].to_string(index=False)); print(topd.to_string(index=False))
