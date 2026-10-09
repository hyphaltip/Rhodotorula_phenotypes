#!/usr/bin/env python3
"""Clonal lineages of the panel, run-adjusted phenotypes, lineage ICC, and GEMMA covariate files.

The 126-strain panel (<= 5 SNP de-cloning, D-49) still consists of a small number of clonal lineages. Lineages here are single-linkage groups of panel
strains with <= LINEAGE_DIFFS pairwise SNP differences (popgen rmuc_core pairwise table). The number of groups is 14 for any cut between 2,000 and
20,000 differences, so 2,000 is taken as the edge of the plateau (decision D-61)."""
from pathlib import Path
import numpy as np, pandas as pd
import warnings; warnings.filterwarnings("ignore")
import statsmodels.formula.api as smf
from scipy import stats
R = Path("analysis/gwas_dh4148/results"); G = R / "gwas"; T = Path("analysis/gwas_dh4148/report/tables")
LINEAGE_DIFFS = 2000
log = lambda m: print(m, flush=True)
ps = pd.read_csv(R / "geno/panel_strains.csv", dtype=str); order = [l.strip() for l in open(R / "geno/panel_samples_in_vcf_order.txt")]
ps = ps.set_index("popgen_strain").loc[order].reset_index(); ps["strain_id"] = ps.strain_id.astype(int)
pw = pd.read_csv("data/raw/popgen-callset-metadata/popgen/variant_qc/declone/rmuc_core.pairwise.tsv.gz", sep="\t"); names = set(ps.popgen_strain)
pw = pw[pw.strain_a.isin(names) & pw.strain_b.isin(names)]
par = {n: n for n in names}
def find(x):
    while par[x] != x: par[x] = par[par[x]]; x = par[x]
    return x
for a, b in zip(pw[pw.diffs <= LINEAGE_DIFFS].strain_a, pw[pw.diffs <= LINEAGE_DIFFS].strain_b): par[find(a)] = find(b)
root = {n: find(n) for n in names}; size = pd.Series(root).value_counts()
lin_id = {r: f"L{i + 1:02d}" for i, r in enumerate(size.index)}
ps["lineage"] = ps.popgen_strain.map(root).map(lin_id); ps["lineage_size"] = ps.lineage.map(ps.lineage.value_counts())
# within / between lineage divergence
ps["k"] = ps.popgen_strain
pwl = pw.copy(); pwl["la"] = pwl.strain_a.map(dict(zip(ps.popgen_strain, ps.lineage))); pwl["lb"] = pwl.strain_b.map(dict(zip(ps.popgen_strain, ps.lineage)))
within = pwl[pwl.la == pwl.lb].diffs; between = pwl[pwl.la != pwl.lb].diffs
log(f"lineages: {ps.lineage.nunique()}; sizes {sorted(ps.lineage.value_counts().tolist(), reverse=True)}")
log(f"pairwise SNP differences: within lineage median {within.median():.0f} (max {within.max():.0f}); between lineages min {between.min():.0f}, median {between.median():.0f}")
ps[["strain_id", "sample_name", "popgen_strain", "lineage", "lineage_size"]].to_csv(R / "lineages.csv", index=False)
ps.groupby("lineage").agg(n_strains=("strain_id", "size"), example_strains=("sample_name", lambda x: ", ".join(list(x)[:4]))).reset_index().to_csv(T / "lineage_sizes.csv", index=False)
# lineage count by cut (re-computed)
rows = []
for thr in (5, 50, 500, 2000, 5000, 20000, 50000):
    pp = {n: n for n in names}
    def f(x):
        while pp[x] != x: pp[x] = pp[pp[x]]; x = pp[x]
        return x
    for a, b in zip(pw[pw.diffs <= thr].strain_a, pw[pw.diffs <= thr].strain_b): pp[f(a)] = f(b)
    comp = pd.Series({n: f(n) for n in names}).value_counts(); rows.append(dict(max_snp_differences_within_group=thr, n_groups=len(comp), largest_groups=", ".join(map(str, comp.values[:6]))))
pd.DataFrame(rows).to_csv(T / "lineages_by_threshold.csv", index=False); log(pd.DataFrame(rows).to_string(index=False))
# runs per strain per metal
w = pd.read_csv("analysis/carotenoid_stress_vs_baseline/results/wells.csv", usecols=["Metal", "run_number", "strain_id"]).drop_duplicates(); w["strain_id"] = w.strain_id.astype(int)
run = {m: w[w.Metal == m].drop_duplicates("strain_id").set_index("strain_id").run_number for m in ("Chromium", "Copper", "Lead")}
ct = []
for m, r in run.items():
    x = ps.assign(run=ps.strain_id.map(r)).dropna(subset=["run"]); tab = pd.crosstab(x.lineage, x.run); tab.insert(0, "Metal", m); ct.append(tab.reset_index())
    chi = stats.chi2_contingency(tab.drop(columns="Metal").values); log(f"{m}: lineage x run chi-square p = {chi[1]:.3f} (n = {len(x)})")
pd.concat(ct).to_csv(T / "lineage_by_run.csv", index=False)
# run-adjusted traits
P = pd.read_csv(G / "pheno_raw.csv", index_col=0); P = P.loc[ps.strain_id.values]
AB = {"Cr": "Chromium", "Cu": "Copper", "Pb": "Lead"}; Padj = P.copy(); icc = []
lin = ps.lineage.values
for t in P.columns:
    y = P[t].copy(); adj = y.copy()
    if t[:2] in AB:
        rr = pd.Series(ps.strain_id.map(run[AB[t[:2]]]).values, index=P.index)
        df = pd.DataFrame(dict(y=y.values, run=rr.values, lin=lin)).dropna()
        m = smf.mixedlm("y ~ C(run)", df, groups=df["lin"]).fit(reml=True, method="lbfgs")
        eff = {r: 0.0 for r in df.run.unique()}
        for k, v in m.fe_params.items():
            if k.startswith("C(run)"): eff[k.split("[T.")[1][:-1]] = v
        mean_eff = np.mean([eff[r] for r in df.run]); adj = y - rr.map(eff).values + mean_eff
    Padj[t] = adj
    for lab, v in (("raw", y), ("run-adjusted", adj)):
        df = pd.DataFrame(dict(y=v.values, lin=lin)).dropna()
        k = df.lin.nunique(); gs = [g.y.values for _, g in df.groupby("lin")]; N = len(df); ni = np.array([len(g) for g in gs])
        grand = df.y.mean(); ssb = sum(n_ * (g.mean() - grand) ** 2 for n_, g in zip(ni, gs)); ssw = sum(((g - g.mean()) ** 2).sum() for g in gs)
        msb = ssb / (k - 1); msw = ssw / (N - k); n0 = (N - (ni ** 2).sum() / N) / (k - 1)
        icc1 = (msb - msw) / (msb + (n0 - 1) * msw); F = msb / msw; p = stats.f.sf(F, k - 1, N - k)
        icc.append(dict(trait=t, version=lab, n_strains=N, n_lineages=k, icc_lineage=icc1, anova_F=F, anova_p=p))
Padj.to_csv(G / "pheno_adj_raw.csv"); icc = pd.DataFrame(icc); icc.to_csv(T / "lineage_icc.csv", index=False)
def rint(x): x = pd.Series(x); return stats.norm.ppf((x.rank() - 0.5) / x.notna().sum())
Z = Padj.apply(rint); Z.to_csv(G / "pheno_rint_adj.csv"); Z.to_csv(G / "pheno_rint_adj.txt", sep=" ", header=False, index=False, na_rep="NA")
# covariates: intercept + lineage dummies (all but the largest lineage)
levels = [l for l in sorted(ps.lineage.unique())][1:]
C = pd.DataFrame({"intercept": 1.0}, index=ps.index)
for l in levels: C[l] = (ps.lineage == l).astype(float).values
C.to_csv(G / "covar_lineage.txt", sep=" ", header=False, index=False)
log(icc[icc.version == "run-adjusted"].round(3).to_string(index=False)); log(icc.pivot(index="trait", columns="version", values="icc_lineage").round(2).to_string())
