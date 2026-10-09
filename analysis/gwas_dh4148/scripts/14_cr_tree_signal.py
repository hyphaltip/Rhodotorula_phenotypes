#!/usr/bin/env python3
"""Cr (and Cu, Pb) phylogenetic signal inside pure haploid R. mucilaginosa: raw against run-adjusted baseline a*, and how run, lineage and tree relate.
Reuses the tree code of the a* report (outgroup-rooted PHYling tree, Pagel's lambda)."""
import sys, warnings; warnings.filterwarnings("ignore")
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, "analysis/carotenoid_stress_vs_baseline/scripts/strat")
src = open("analysis/carotenoid_stress_vs_baseline/scripts/strat/s2_phylo.py").read().split("w = load_wells(); rows = []")[0]
ns = {"__name__": "s2lib"}; exec(compile(src, "s2_phylo_head", "exec"), ns)
C, name2i, pagel = ns["C"], ns["name2i"], ns["pagel"]
from strat_common import *
T = Path("analysis/gwas_dh4148/report/tables")
w = load_wells(); w["strain_id"] = w.strain_id.astype(int)
lin = pd.read_csv("analysis/gwas_dh4148/results/lineages.csv")[["strain_id", "lineage"]]
rows = []; detail = []
for m in ["Chromium", "Copper", "Lead"]:
    g = w[(w.Metal == m) & (w.species == "Rhodotorula mucilaginosa") & (w.conc == 0)]
    b = g.groupby("strain_id").agg(a=("a", "mean"), run=("run_number", "first"), tip=("tip", "first")).reset_index().merge(lin, on="strain_id", how="left")
    b = b[b.tip.notna() & b.tip.isin(name2i)].copy()
    b["a_runadj"] = b.a - b.groupby("run").a.transform("mean") + b.a.mean()
    ix = np.array([name2i[t] for t in b.tip]); Cs = C[np.ix_(ix, ix)] + 1e-6 * np.mean(np.diag(C)) * np.eye(len(ix))
    for lab, col in (("baseline a*", "a"), ("baseline a*, run effect removed", "a_runadj")):
        lam, ll, p0, p1 = pagel(b[col].values, Cs); rows.append(dict(Metal=m, trait=lab, n_strains=len(b), n_runs=b.run.nunique(), pagel_lambda=lam, p_lambda_gt_0=p0))
    # variance of baseline a* explained by run, by lineage
    for fac in ("run", "lineage"):
        x = b.dropna(subset=[fac]); k = x[fac].nunique(); grand = x.a.mean(); ssb = sum(len(v) * (v.a.mean() - grand) ** 2 for _, v in x.groupby(fac)); sst = ((x.a - grand) ** 2).sum()
        detail.append(dict(Metal=m, factor=fac, n_groups=k, n_strains=len(x), r2=ssb / sst, adj_r2=1 - (1 - ssb / sst) * (len(x) - 1) / (len(x) - k)))
    if m == "Chromium":
        mm = b.groupby("run").agg(n=("a", "size"), mean_a=("a", "mean"), sd_a=("a", "std")).round(2); print(mm)
        print(b.groupby("lineage").agg(n=("a", "size"), mean_a=("a", "mean"), n_runs=("run", "nunique")).round(2).sort_values("n", ascending=False).head(8))
        print(pd.crosstab(b.lineage, b.run).loc[lambda d: d.sum(1) >= 5])
r = pd.DataFrame(rows); r.to_csv(T / "cr_tree_signal_run_adjusted.csv", index=False); d = pd.DataFrame(detail); d.to_csv(T / "baseline_a_variance_by_run_and_lineage.csv", index=False)
print(r.round(3).to_string(index=False)); print(d.round(3).to_string(index=False))
