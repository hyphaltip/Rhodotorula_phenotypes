#!/usr/bin/env python3
"""Check the phylogenetic-signal numbers: dependence on the diagonal jitter, and on collapsing near-identical tips."""
import sys, re
import numpy as np, pandas as pd
src = open("analysis/carotenoid_stress_vs_baseline/scripts/strat/s2_phylo.py").read().split("w = load_wells(); rows = []")[0]
ns = {"__name__": "s2lib"}; exec(compile(src, "s2_phylo_head", "exec"), ns)
for k in ("C", "name2i", "pagel", "K_perm", "n", "depth", "tips", "tip_ix", "par", "ln", "kids"): globals()[k] = ns[k]
sys.path.insert(0, "analysis/carotenoid_stress_vs_baseline/scripts/strat")
from strat_common import *
term = np.array([ln[v] for v in tips]); print("terminal branch lengths: min %.2e, p10 %.2e, median %.2e, max %.2e; zero-length tips: %d of %d" % (term.min(), np.quantile(term, .1), np.median(term), term.max(), (term == 0).sum(), len(term)))
d = np.add.outer(np.diag(C), np.diag(C)) - 2 * C; np.fill_diagonal(d, np.inf); nn = d.min(1)
print("nearest-neighbour patristic distance: min %.2e, p10 %.2e, median %.2e; tips with a neighbour closer than 1e-5: %d" % (nn.min(), np.quantile(nn, .1), np.median(nn), (nn < 1e-5).sum()))
print("mean diag(C) = %.4f" % np.mean(np.diag(C)))
w = load_wells(); rng = np.random.default_rng(3); rows = []
for m in ("Copper", "Lead", "Chromium"):
    g = w[w.Metal == m].copy(); g["a_adj"], _ = size_adjusted(g)
    b = g[g.conc == 0].groupby("strain_id").agg(a_adj=("a_adj", "mean"), tip=("tip", "first")); s = b[b.tip.notna() & b.tip.isin(name2i)]
    ix = np.array([name2i[t] for t in s.tip]); y = s.a_adj.values; Cs0 = C[np.ix_(ix, ix)]
    for jit in (1e-8, 1e-6, 1e-4, 1e-3, 1e-2):
        Cs = Cs0 + jit * np.mean(np.diag(C)) * np.eye(len(ix))
        lam, ll, p0, p1 = pagel(y, Cs); K, pK = K_perm(y, Cs, B=199)
        rows.append(dict(Metal=m, jitter=jit, n=len(y), cond=np.linalg.cond(Cs), lam=lam, p_lam0=p0, K=K, pK=pK)); print(rows[-1], flush=True)
    # non-tree check: Moran-like -- does phenotype similarity track patristic distance? Spearman(patristic distance, |dy|) over all pairs
    D = np.add.outer(np.diag(Cs0), np.diag(Cs0)) - 2 * Cs0; iu = np.triu_indices(len(ix), 1)
    from scipy import stats as st
    print(m, "Spearman(patristic distance, |difference in a*adj|) over pairs:", st.spearmanr(D[iu], np.abs(y[:, None] - y[None, :])[iu])[0])
pd.DataFrame(rows).to_csv(REP / "tables/s2_jitter_sensitivity.csv", index=False)
