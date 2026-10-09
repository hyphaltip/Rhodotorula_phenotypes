#!/usr/bin/env python3
"""Summarise the GEMMA scans: inflation with and without the relatedness matrix, Manhattan and QQ plots, lead SNPs per locus, annotation."""
import gzip, sys
sys.path.insert(0, "analysis/gwas_dh4148/scripts")
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy import stats

from geno_lib import *

G = R / "gwas"; T = Path("analysis/gwas_dh4148/report/tables"); F = Path("analysis/gwas_dh4148/report/figures")
OUTD = sys.argv[1] if len(sys.argv) > 1 else "output"; PHENO = sys.argv[2] if len(sys.argv) > 2 else "pheno_rint.csv"; TAG = sys.argv[3] if len(sys.argv) > 3 else ""   # TAG is appended to every file name
LINCOV = len(sys.argv) > 4 and sys.argv[4] == "lineage"
MEFF = sum(1 for _ in open(GENO / "pruned_r2_0.5.snps.txt")); BONF = 0.05 / MEFF; SUGG = 1.0 / MEFF   # may be replaced below by unique patterns
log = lambda m: print(m, flush=True)
tl = pd.read_csv(G / "trait_list.csv"); traits = tl.trait.tolist()
ids, X, sid, names = load_geno(); n = X.shape[1]
# effective number of tests = number of distinct genotype patterns (perfectly linked SNPs collapse). With lineage covariates the pattern is the
# within-lineage residual, and SNPs that do not vary inside any lineage are not testable.
if LINCOV:
    lin = pd.read_csv(R / "lineages.csv").set_index("strain_id").loc[sid, "lineage"].values; Xr = X.copy()
    for l in np.unique(lin): Xr[:, lin == l] = X[:, lin == l] - np.nanmean(X[:, lin == l], axis=1, keepdims=True)
    Xr = np.round(np.nan_to_num(Xr, nan=9), 2); var_in = (np.abs(Xr) > 1e-6).any(1) & (Xr != 9).any(1); pat = Xr[var_in]
else: pat = np.where(np.isnan(X), -1, X).astype(np.int8)
MEFF = len(np.unique(pat, axis=0)); BONF = 0.05 / MEFF; SUGG = 1.0 / MEFF
log(f"distinct genotype patterns (effective number of tests): {MEFF:,} of {len(ids):,} SNPs; r2 < 0.5 pruned set: {sum(1 for _ in open(GENO / 'pruned_r2_0.5.snps.txt'))}")
chrom, pos = split_id(ids); ct = contigs(); off = dict(zip(ct.chr, np.concatenate([[0], np.cumsum(ct.len.values)[:-1]] ))); keep = np.isin(chrom, ct.chr)
log(f"Bonferroni p < {BONF:.2e}; suggestive p < {SUGG:.2e}")
P = pd.read_csv(G / PHENO, index_col=0)
# naive OLS scan (no relatedness) for the inflation comparison
Xm = np.where(np.isnan(X), np.nanmean(X, axis=1, keepdims=True), X); Xc = Xm - Xm.mean(1, keepdims=True); sx = (Xc ** 2).sum(1)
def lam(p): p = np.asarray(p); p = p[np.isfinite(p) & (p > 0)]; return np.median(stats.chi2.isf(p, 1)) / 0.4549
res = {}; rows = []; lead_rows = []
for t in traits:
    a = pd.read_csv(G / OUTD / f"{t}.assoc.txt.gz", sep="\t"); a["chr"], a["pos"] = split_id(a.rs.values); res[t] = a
    y = P[t].values; ok = ~np.isnan(y); yc = y[ok] - y[ok].mean()
    b = (Xc[:, ok] @ yc) / sx; rss = ((yc ** 2).sum() - b ** 2 * sx) ; se = np.sqrt(rss / (ok.sum() - 2) / sx); pn = 2 * stats.t.sf(np.abs(b / se), ok.sum() - 2)
    pv = float(tl.set_index("trait").loc[t, "n_strains"])
    log_t = (G / OUTD / f"{t}.log.txt").read_text(); pve = float([l for l in log_t.splitlines() if "pve estimate" in l][0].split("=")[1]); pvese = float([l for l in log_t.splitlines() if "se(pve)" in l][0].split("=")[1])
    rows.append(dict(trait=t, n_strains=int(pv), n_snps=len(a), lambda_naive_ols=lam(pn), lambda_lmm=lam(a.p_wald), snp_pve=pve, snp_pve_se=pvese, min_p_lmm=a.p_wald.min(), n_bonferroni=int((a.p_wald < BONF).sum()), n_suggestive=int((a.p_wald < SUGG).sum()), n_suggestive_expected_null=len(a) * SUGG))
inf = pd.DataFrame(rows); inf.to_csv(T / f"gwas_inflation_pve{TAG}.csv", index=False)
# loci: clump suggestive SNPs within 50 kb, lead = min p
genes = load_genes()
def near_gene(c, p):
    g = genes[genes.chr == c]; inside = g[(g.start <= p) & (g.end >= p)]
    if len(inside): r = inside.iloc[0]; return r.gene_id, r.symbol, r["product"], 0
    d = np.minimum(abs(g.start - p), abs(g.end - p)); 
    if not len(g): return "", "", "", np.nan
    i = d.values.argmin(); r = g.iloc[i]; return r.gene_id, r.symbol, r["product"], int(d.values[i])
for t, a in res.items():
    s = a[a.p_wald < SUGG].sort_values("p_wald")
    used = []
    for _, r in s.iterrows():
        if any(r.chr == c and abs(r.pos - p) < 50000 for c, p in used): continue
        used.append((r.chr, r.pos)); cl = s[(s.chr == r.chr) & (abs(s.pos - r.pos) < 50000)]
        lead_rows.append(dict(trait=t, lead_snp=r.rs, chr=r.chr, pos=r.pos, af=r.af, beta=r.beta, se=r.se, p_wald=r.p_wald, bonferroni_sig=bool(r.p_wald < BONF), n_suggestive_snps_in_locus=len(cl), locus_start=cl.pos.min(), locus_end=cl.pos.max()))
lead = pd.DataFrame(lead_rows)
if len(lead):
    an = ann_for(lead.lead_snp); lead[["effect", "impact", "ann_gene", "hgvs_p"]] = [an.get(k, ("", "", "", "")) for k in lead.lead_snp]
    ng = [near_gene(c, p) for c, p in zip(lead.chr, lead.pos)]; lead[["gene_id", "gene_symbol", "product", "distance_bp"]] = pd.DataFrame(ng, index=lead.index)
    lead["n_genes_in_locus"] = [int(((genes.chr == c) & (genes.end >= s0) & (genes.start <= e0)).sum()) for c, s0, e0 in zip(lead.chr, lead.locus_start, lead.locus_end)]
lead.sort_values("p_wald").to_csv(T / f"gwas_loci{TAG}.csv", index=False)
top = pd.concat([a.sort_values("p_wald").head(10).assign(trait=t) for t, a in res.items()])[["trait", "rs", "chr", "pos", "af", "beta", "se", "p_wald"]]; top.to_csv(T / f"gwas_top10_per_trait{TAG}.csv", index=False)
log(inf.round(3).to_string(index=False)); log(f"loci at p < {SUGG:.1e}: {len(lead)}; Bonferroni-significant loci: {int(lead.bonferroni_sig.sum()) if len(lead) else 0}")
# figures: Manhattan per group of traits, QQ grid
def xpos(a): return a.pos.values + a.chr.map(off).fillna(0).values
groups = {"dose0": [t for t in traits if t.startswith("d0")], "Cr": [t for t in traits if t.startswith("Cr")], "Cu": [t for t in traits if t.startswith("Cu")], "Pb": [t for t in traits if t.startswith("Pb")]}
for gname, ts in groups.items():
    fig, ax = plt.subplots(len(ts), 1, figsize=(13, 2.2 * len(ts) + 0.6), sharex=True, squeeze=False)
    for k, t in enumerate(ts):
        a = res[t]; a = a[a.chr.isin(ct.chr)]; col = np.where(a.chr.map({c: i for i, c in enumerate(ct.chr)}) % 2 == 0, "#0072B2", "#56B4E9")
        ax[k, 0].scatter(xpos(a), -np.log10(a.p_wald), s=2, c=col, rasterized=True); ax[k, 0].axhline(-np.log10(BONF), color="r", lw=.7); ax[k, 0].axhline(-np.log10(SUGG), color="grey", ls=":", lw=.7)
        ax[k, 0].set_ylabel(t.replace("_", " "), fontsize=7); ax[k, 0].set_ylim(0, max(5, -np.log10(a.p_wald.min()) + .3))
    ax[-1, 0].set_xlabel("position along concatenated contigs >= 100 kb (red: Bonferroni over distinct genotype patterns; dotted: 1 expected false positive)")
    fig.tight_layout(); fig.savefig(F / f"gwas_manhattan_{gname}{TAG}.png", dpi=150); plt.close(fig)
fig, ax = plt.subplots(4, 5, figsize=(15, 11))
for k, t in enumerate(traits):
    a = ax.flat[k]; p = np.sort(res[t].p_wald.values); e = (np.arange(1, len(p) + 1) - .5) / len(p)
    a.scatter(-np.log10(e), -np.log10(p), s=2, rasterized=True); m = max(-np.log10(e)); a.plot([0, m], [0, m], "r", lw=.7); a.set_title(f"{t} (lambda {inf.set_index('trait').loc[t, 'lambda_lmm']:.2f})", fontsize=7)
for a in ax.flat[len(traits):]: a.axis("off")
fig.tight_layout(); fig.savefig(F / f"gwas_qq{TAG}.png", dpi=150); plt.close(fig)
fig, ax = plt.subplots(figsize=(9, 4.4)); x = np.arange(len(inf)); ax.bar(x - .2, inf.lambda_naive_ols, .4, label="OLS, no relatedness"); ax.bar(x + .2, inf.lambda_lmm, .4, label="LMM with relatedness matrix"); ax.axhline(1, color="k", lw=.7)
ax.set_xticks(x); ax.set_xticklabels(inf.trait, rotation=70, fontsize=7); ax.set_ylabel("genomic inflation lambda_GC"); ax.legend(fontsize=8); fig.tight_layout(); fig.savefig(F / f"gwas_inflation{TAG}.png", dpi=160); plt.close(fig)
# relatedness structure of the panel
K = np.loadtxt(G / "output" / "K.cXX.txt"); w, v = np.linalg.eigh(K); o = np.argsort(w)[::-1]; w, v = w[o], v[:, o]
sp = pd.read_csv("data/metadata/strain-curation/strain_curation.csv", dtype=str).set_index("strain_id")
fig, ax = plt.subplots(1, 2, figsize=(11, 4.4)); ax[0].scatter(v[:, 0] * np.sqrt(w[0]), v[:, 1] * np.sqrt(w[1]), s=14); ax[0].set_xlabel(f"PC1 ({100 * w[0] / w.sum():.1f}% of relatedness)"); ax[0].set_ylabel(f"PC2 ({100 * w[1] / w.sum():.1f}%)")
ax[1].bar(range(1, 21), 100 * w[:20] / w.sum()); ax[1].set_xlabel("PC"); ax[1].set_ylabel("% of relatedness variance"); fig.suptitle("Genetic structure of the 126-strain panel (centred relatedness matrix, 244,956 SNPs)", fontsize=9)
fig.tight_layout(); fig.savefig(F / "gwas_panel_structure.png", dpi=160); plt.close(fig)
pd.DataFrame({"pc": range(1, 11), "pct_var": 100 * w[:10] / w.sum()}).to_csv(T / "gwas_panel_pc_variance.csv", index=False)
np.save(G / "K_eig_vals.npy", w)
