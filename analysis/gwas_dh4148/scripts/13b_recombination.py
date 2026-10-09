#!/usr/bin/env python3
"""Is there recombination inside the large clonal lineages (L01-L03)? Haploid genotypes of informative SNPs (minor allele count >= 2).
 (1) Four-gamete test and r2 for SNP pairs by physical distance, and for pairs on different contigs. Under strict clonality every pair is equally linked, so
     the incompatible fraction does not depend on distance. Under recombination, close pairs stay compatible and distant or unlinked pairs become incompatible.
 (2) Windowed genetic distances: correlation of strain-by-strain distance matrices between 300 kb windows (adjacent windows on one contig, distant windows, other contigs)."""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy import stats
R = Path("analysis/gwas_dh4148/results/recomb"); T = Path("analysis/gwas_dh4148/report/tables"); F = Path("analysis/gwas_dh4148/report/figures")
rng = np.random.default_rng(11); log = lambda m: print(m, flush=True)
BINS = [(0, 1e3), (1e3, 5e3), (5e3, 2e4), (2e4, 1e5), (1e5, 5e5), (5e5, 1e9)]; LAB = ["<1 kb", "1-5 kb", "5-20 kb", "20-100 kb", "100-500 kb", ">500 kb (same contig)"]
rows = []; wrows = []
fig, ax = plt.subplots(1, 3, figsize=(15, 4.4), sharey=True)
for li, L in enumerate(["L01", "L02", "L03"]):
    d = pd.read_csv(R / f"{L}.gt.tsv.gz", sep="\t", header=None, na_values="NA"); chrom = d[0].values; pos = d[1].values; G = d.iloc[:, 4:].values.astype(float).T   # n x M
    n, M = G.shape; log(f"{L}: {n} strains, {M} informative SNPs, {len(set(chrom))} contigs")
    if M < 50: continue
    maf = np.nanmean(G, axis=0); G[:, maf > 0.5] = 1 - G[:, maf > 0.5]   # minor allele = 1
    # pair sampling: stratified by distance class; "other contig" = random pairs on different contigs
    cls = {}
    idx_by_c = {c: np.where(chrom == c)[0] for c in np.unique(chrom)}
    def sample_same(lo, hi, k):
        out = []; tries = 0
        while len(out) < k and tries < 60:
            a = rng.integers(0, M, 20000); c = chrom[a]; b = np.empty_like(a)
            for cc in np.unique(c):
                m = c == cc; ia = idx_by_c[cc]; b[m] = rng.choice(ia, m.sum())
            dist = np.abs(pos[a] - pos[b]); ok = (dist >= lo) & (dist < hi) & (a != b); out += list(zip(a[ok], b[ok])); tries += 1
        return np.array(out[:k])
    K = 30000
    for (lo, hi), lab in zip(BINS, LAB): cls[lab] = sample_same(lo, hi, K)
    a = rng.integers(0, M, 200000); b = rng.integers(0, M, 200000); ok = chrom[a] != chrom[b]; cls["other contig"] = np.column_stack([a[ok], b[ok]])[:K]
    for lab, pr in cls.items():
        if len(pr) == 0: continue
        x = G[:, pr[:, 0]]; y = G[:, pr[:, 1]]; v = ~np.isnan(x) & ~np.isnan(y); nn = v.sum(0)
        x = np.where(v, x, 0); y = np.where(v, y, 0)
        n11 = (x * y).sum(0); n10 = (x * (1 - y) * v).sum(0); n01 = ((1 - x) * y * v).sum(0); n00 = ((1 - x) * (1 - y) * v).sum(0)
        px = (n11 + n10) / nn; py = (n11 + n01) / nn; D = n11 / nn - px * py; den = px * (1 - px) * py * (1 - py); r2 = np.where(den > 0, D * D / den, np.nan)
        inc1 = (np.minimum.reduce([n11, n10, n01, n00]) >= 1); inc2 = (np.minimum.reduce([n11, n10, n01, n00]) >= 2)
        rows.append(dict(lineage=L, n_strains=n, n_snps=M, pair_class=lab, n_pairs=len(pr), mean_r2=np.nanmean(r2), frac_four_gametes=inc1.mean(), frac_four_gametes_each_ge2=inc2.mean()))
    # windowed distance-matrix correlations
    W = 300000; wins = {}
    for c in np.unique(chrom):
        ii = idx_by_c[c]; wb = pos[ii] // W
        for w_ in np.unique(wb):
            s = ii[wb == w_]
            if len(s) >= 15: wins[(c, int(w_))] = s
    def dm(s):
        g = G[:, s]; nan = np.isnan(g); g0 = np.where(nan, 0, g); cnt = (~nan).astype(float)
        diff = (g0[:, None, :] != g0[None, :, :]) & (cnt[:, None, :] * cnt[None, :, :] > 0); return diff.sum(2)[np.triu_indices(n, 1)].astype(float)
    keys = list(wins)
    if len(keys) < 3: continue
    D = {k: dm(wins[k]) for k in keys}; res = {"adjacent windows (same contig)": [], "distant windows (same contig, > 1 Mb at 300 kb windows)": [], "different contigs": []}
    for _ in range(4000):
        k1, k2 = rng.choice(len(keys), 2, replace=False); a_, b_ = keys[k1], keys[k2]
        if D[a_].std() == 0 or D[b_].std() == 0: continue
        rho = stats.spearmanr(D[a_], D[b_])[0]
        if a_[0] != b_[0]: res["different contigs"].append(rho)
        elif abs(a_[1] - b_[1]) == 1: res["adjacent windows (same contig)"].append(rho)
        elif abs(a_[1] - b_[1]) > 3: res["distant windows (same contig, > 1 Mb at 300 kb windows)"].append(rho)
    for k, v in res.items(): wrows.append(dict(lineage=L, window_pair_class=k, n_window_pairs=len(v), mean_spearman=np.mean(v) if v else np.nan, q05=np.quantile(v, .05) if v else np.nan, n_windows=len(keys)))
    s = pd.DataFrame([r for r in rows if r["lineage"] == L]); ax[li].plot(range(len(s)), s.frac_four_gametes, "-o", label="all four gametes (>=1)"); ax[li].plot(range(len(s)), s.frac_four_gametes_each_ge2, "-s", label="each gamete >= 2 strains")
    ax[li].set_xticks(range(len(s))); ax[li].set_xticklabels(s.pair_class, rotation=40, ha="right", fontsize=7); ax[li].set_title(f"{L}: {n} strains, {M:,} informative SNPs", fontsize=9)
ax[0].set_ylabel("fraction of SNP pairs with all four gametes"); ax[0].legend(fontsize=7); fig.tight_layout(); fig.savefig(F / "recombination_four_gamete.png", dpi=160); plt.close(fig)
pd.DataFrame(rows).to_csv(T / "recombination_pair_classes.csv", index=False); pd.DataFrame(wrows).to_csv(T / "recombination_window_distance_correlation.csv", index=False)
log(pd.DataFrame(rows).round(3).to_string(index=False)); log(pd.DataFrame(wrows).round(3).to_string(index=False))
