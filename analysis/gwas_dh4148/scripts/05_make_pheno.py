#!/usr/bin/env python3
"""Phenotype matrix for the 126-strain GWAS panel (rows in the genotype order). Rank-based inverse normal transform per trait.
Writes results/gwas/pheno_raw.csv, pheno_rint.txt (GEMMA, no header) and trait_list.csv."""
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
R = Path("analysis/gwas_dh4148/results"); G = R / "gwas"; G.mkdir(exist_ok=True, parents=True)
ps = pd.read_csv(R / "geno/panel_strains.csv", dtype=str); order = [l.strip() for l in open(R / "geno/panel_samples_in_vcf_order.txt")]
assert set(order) == set(ps.popgen_strain) and len(order) == 126
ps = ps.set_index("popgen_strain").loc[order].reset_index(); sid = ps.strain_id.astype(int).values
tr = pd.read_csv(R / "strain_traits.csv")
bl = pd.read_csv("analysis/gwas_dh4148/report/tables/dose0_strain_blups.csv").pivot(index="strain_id", columns="trait", values="blup")
P = pd.DataFrame(index=sid)
for tname, c in (("d0_a", "a"), ("d0_lnA", "lnA"), ("d0_rgr", "rgr")): P[tname] = bl[c].reindex(sid).values
M = {"Chromium": "Cr", "Copper": "Cu", "Lead": "Pb"}
for m, ab in M.items():
    t = tr[tr.Metal == m].set_index("strain_id").reindex(sid)
    P[f"{ab}_da_top"] = t.da_top.values; P[f"{ab}_a_auc"] = t.a_auc.values; P[f"{ab}_relarea_top"] = t.relarea_top.values
    P[f"{ab}_relarea_auc"] = t.relarea_auc.values; P[f"{ab}_relrgr_auc"] = t.relrgr_auc.values
    if m != "Copper": P[f"{ab}_logIC50"] = np.log(t.ic50.values)      # Copper: 46% of panel strains are right-censored (relative area >= 0.5), IC50 not used
P = P.dropna(axis=1, how="all"); P.index.name = "strain_id"
def rint(x):
    x = pd.Series(x); r = x.rank(); n = x.notna().sum(); return stats.norm.ppf((r - 0.5) / n)
Z = P.apply(rint)
P.to_csv(G / "pheno_raw.csv"); Z.to_csv(G / "pheno_rint.csv")
Z.to_csv(G / "pheno_rint.txt", sep=" ", header=False, index=False, na_rep="NA")
tl = pd.DataFrame({"idx": range(1, P.shape[1] + 1), "trait": P.columns, "n_strains": P.notna().sum().values, "mean_raw": P.mean().values, "sd_raw": P.std().values})
tl.to_csv(G / "trait_list.csv", index=False); print(tl.round(3).to_string(index=False))
