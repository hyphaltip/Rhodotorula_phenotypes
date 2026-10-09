#!/usr/bin/env python3
"""Clone-corrected association inside each large lineage (L01, L02, L03), separately.
Within a lineage all strains share one background (<= 2,010 SNP differences), so the informative SNPs are lineage-private mutations (minor allele count >= 2).
Per trait: OLS of the run-adjusted rank-normalised trait on each SNP; family-wise p from 2,000 permutations of the trait across strains (max |t| over SNPs)."""
import sys
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, "analysis/gwas_dh4148/scripts")
from geno_lib import best_ann
R = Path("analysis/gwas_dh4148/results"); T = Path("analysis/gwas_dh4148/report/tables"); rng = np.random.default_rng(5); NPERM = 2000
log = lambda m: print(m, flush=True)
lin = pd.read_csv(R / "lineages.csv"); Z = pd.read_csv(R / "gwas/pheno_rint_adj.csv", index_col=0)
out = []; summ = []
for L in ("L01", "L02", "L03"):
    samples = [l.strip() for l in open(R / f"recomb/{L}.samples.txt")]; names = [s[: len(s) // 2] for s in samples]
    sid = lin.set_index("popgen_strain").loc[names, "strain_id"].values
    d = pd.read_csv(R / f"recomb/{L}.gt.tsv.gz", sep="\t", header=None, na_values="NA"); an = pd.read_csv(R / f"recomb/{L}.ann.tsv.gz", sep="\t", header=None, names=["c", "p", "r", "a", "ann"], dtype=str)
    X = d.iloc[:, 4:].values.astype(float).T; chrom = d[0].values; pos = d[1].values; ref = d[2].values; alt = d[3].values; ann = an.ann.fillna("").values
    for t in Z.columns:
        y = Z.loc[sid, t].values; ok = ~np.isnan(y)
        if ok.sum() < 12: continue
        Xo = X[ok]; yo = y[ok]; n = ok.sum(); Xm = np.where(np.isnan(Xo), np.nanmean(Xo, 0, keepdims=True), Xo); Xc = Xm - Xm.mean(0); sx = (Xc ** 2).sum(0); good = sx > 1e-9
        def tstat(yv):
            yc = yv - yv.mean(); b = (Xc.T @ yc) / np.where(good, sx, 1); rss = (yc ** 2).sum() - b * b * sx; se = np.sqrt(np.maximum(rss, 1e-12) / (n - 2) / np.where(good, sx, 1)); return np.where(good, b / se, 0), b
        tv, b = tstat(yo); from scipy import stats
        p = 2 * stats.t.sf(np.abs(tv), n - 2); mx = np.array([np.abs(tstat(rng.permutation(yo))[0]).max() for _ in range(NPERM)])
        padj = np.array([(mx >= abs(v)).mean() for v in tv]); padj = np.maximum(padj, 1 / (NPERM + 1))
        npat = len({tuple(np.round(r, 3)) for r in Xm.T[good]})
        summ.append(dict(lineage=L, trait=t, n_strains=int(n), n_snps=int(good.sum()), distinct_patterns=npat, min_p=p[good].min(), min_p_family_wise=padj[good].min(), n_family_wise_lt_0_05=int((padj[good] < 0.05).sum())))
        for i in np.argsort(np.where(good, p, 2))[:10]:
            e = best_ann(ann[i], alt[i]) if ann[i] else ("", "", "", "")
            out.append(dict(lineage=L, trait=t, chr=chrom[i], pos=int(pos[i]), ref=ref[i], alt=alt[i], n_alt=int(np.nansum(Xo[:, i])), beta=b[i], p=p[i], p_family_wise=padj[i], effect=e[0], impact=e[1], gene_id=e[2], hgvs_p=e[3]))
    log(f"{L}: {len(d)} SNPs, {len(sid)} strains done")
S = pd.DataFrame(summ); S.to_csv(T / "within_lineage_scan_summary.csv", index=False); O = pd.DataFrame(out); O.to_csv(T / "within_lineage_scan_top_snps.csv", index=False)
log(S[S.trait.str.startswith("Cr")].round(4).to_string(index=False)); log(f"traits x lineages with family-wise p < 0.05: {(S.min_p_family_wise < 0.05).sum()} of {len(S)}")
