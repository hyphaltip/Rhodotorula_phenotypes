#!/usr/bin/env python3
"""Selection on baseline without regression to the mean.
Strains need >= 2 replicate wells at dose 0 and at the top dose. Replicate wells are split at random into halves A and B (1000 splits).
Naive: select tertiles on the baseline from all wells, measure change (top dose - 0) from the same wells.
Split: select tertiles on the baseline in half A; measure the change in half B.
Also: reliability of strain-specific change = Spearman(change in A, change in B)."""
import sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, "analysis/carotenoid_stress_vs_baseline/scripts/strat")
from strat_common import *
rng = np.random.default_rng(11); NSPLIT = 1000; NBOOT = 300
w = load_wells(species=False); rows = []; fig, ax = plt.subplots(2, 4, figsize=(19, 8.4)); store = {}
for k, m in enumerate(["Chromium", "Copper", "Iron", "Lead"]):
    g = w[w.Metal == m].copy(); g["a_adj"], _ = size_adjusted(g); dm = top_dose(g)
    for meas in ("a", "a_adj"):
        cnt = g[g.conc.isin([0, dm])].groupby(["strain_id", "conc"]).size().unstack()
        ok = cnt[(cnt[0] >= 2) & (cnt[dm] >= 2)].index
        v = g[g.strain_id.isin(ok) & g.conc.isin([0, dm])][["strain_id", "conc", meas]].values
        d = {(s, c): [] for s in ok for c in (0, dm)}
        for s, c, x in v: d[(s, c)].append(x)
        S = list(ok); nS = len(S)
        def gap(base, chg):
            lo, hi = np.quantile(base, [1/3, 2/3]); return chg[base >= hi].mean() - chg[base <= lo].mean()
        b_all = np.array([np.mean(d[(s, 0)]) for s in S]); c_all = np.array([np.mean(d[(s, dm)]) for s in S]) - b_all
        naive = gap(b_all, c_all); rho_naive = stats.spearmanr(b_all, c_all)[0]
        def one(idx):
            A0, B0, AT, BT = [], [], [], []
            for s in idx:
                a0 = rng.permutation(d[(s, 0)]); aT = rng.permutation(d[(s, dm)]); h0 = len(a0) // 2; hT = len(aT) // 2
                A0.append(a0[:h0].mean()); B0.append(a0[h0:].mean()); AT.append(aT[:hT].mean()); BT.append(aT[hT:].mean())
            A0, B0, AT, BT = map(np.array, (A0, B0, AT, BT)); dA, dB = AT - A0, BT - B0
            return gap(A0, dB), stats.spearmanr(dA, dB)[0], stats.spearmanr(A0, B0)[0], stats.spearmanr(A0, dB)[0], (A0, dB)
        sp = [one(S) for _ in range(NSPLIT)]
        sp_gap = np.array([x[0] for x in sp]); rel = np.array([x[1] for x in sp]); rb = np.array([x[2] for x in sp]); rAB = np.array([x[3] for x in sp])
        # strain bootstrap for CIs of naive gap and one split gap
        bn, bs = [], []
        for _ in range(NBOOT):
            idx = rng.integers(0, nS, nS); bn.append(gap(b_all[idx], c_all[idx]))
            bs.append(one([S[i] for i in idx])[0])  # strains repeated in the bootstrap sample share wells; splits re-drawn
        rows.append(dict(Metal=m, measure=meas, top_dose=dm, n_strains=nS,
            naive_gap=naive, naive_ci_lo=np.percentile(bn, 2.5), naive_ci_hi=np.percentile(bn, 97.5), naive_spearman_baseline_vs_change=rho_naive,
            split_gap_mean=sp_gap.mean(), split_gap_ci_lo=np.percentile(bs, 2.5), split_gap_ci_hi=np.percentile(bs, 97.5), split_spearman_baselineA_vs_changeB=rAB.mean(),
            reliability_change_A_vs_B=rel.mean(), reliability_ci_lo=np.percentile(rel, 2.5), reliability_ci_hi=np.percentile(rel, 97.5), reliability_baseline_A_vs_B=rb.mean()))
        print(rows[-1], flush=True)
        if meas == "a_adj":
            store[m] = (b_all, c_all, one(S)[4])
tab = pd.DataFrame(rows); tab.to_csv(REP / "tables/s5_split_half.csv", index=False)
for k, m in enumerate(["Chromium", "Copper", "Iron", "Lead"]):
    b, c, (A0, dB) = store[m]
    for r_, (x, y, ttl) in enumerate(((b, c, "naive: baseline and change from the same wells"), (A0, dB, "split-half: baseline from half A, change from half B"))):
        ax[r_, k].scatter(x, y, s=12, alpha=.6, color=COL[m]); sl, ic, *_ = stats.linregress(x, y); xx = np.array([x.min(), x.max()]); ax[r_, k].plot(xx, ic + sl * xx, "k-", lw=1)
        ax[r_, k].set_title(f"{m}: {ttl}\nSpearman {stats.spearmanr(x, y)[0]:.2f}, slope {sl:.2f}", fontsize=8); ax[r_, k].set_xlabel("baseline size-adjusted a* (0 dose)"); ax[r_, k].set_ylabel("change in a* at top dose")
fig.suptitle("Regression to the mean in the baseline-vs-change relationship (size-adjusted a*; one random split shown)", fontsize=10)
fig.tight_layout(); fig.savefig(REP / "figures/s5_split_half_scatter.png", dpi=170); plt.close(fig)
t = tab[tab.measure == "a_adj"]; fig, ax = plt.subplots(1, 2, figsize=(11, 4.2)); x = np.arange(len(t))
ax[0].bar(x - .2, t.naive_gap, .38, yerr=[t.naive_gap - t.naive_ci_lo, t.naive_ci_hi - t.naive_gap], color="#D55E00", label="naive"); ax[0].bar(x + .2, t.split_gap_mean, .38, yerr=[t.split_gap_mean - t.split_gap_ci_lo, t.split_gap_ci_hi - t.split_gap_mean], color="#0072B2", label="split-half")
ax[0].axhline(0, color="k", lw=.6); ax[0].set_xticks(x); ax[0].set_xticklabels(t.Metal); ax[0].set_ylabel("change in a*: top-baseline tertile minus bottom tertile"); ax[0].legend(); ax[0].set_title("Selection gap, size-adjusted a*", fontsize=9)
ax[1].bar(x, t.reliability_change_A_vs_B, yerr=[t.reliability_change_A_vs_B - t.reliability_ci_lo, t.reliability_ci_hi - t.reliability_change_A_vs_B], color="#009E73"); ax[1].axhline(0, color="k", lw=.6)
ax[1].set_xticks(x); ax[1].set_xticklabels(t.Metal); ax[1].set_ylabel("Spearman, change in half A vs half B"); ax[1].set_title("Reliability of the strain-specific change", fontsize=9)
fig.tight_layout(); fig.savefig(REP / "figures/s5_split_half_summary.png", dpi=170); plt.close(fig)
print("done")
