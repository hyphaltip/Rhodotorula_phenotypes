#!/usr/bin/env python3
"""Non-additive genetic architecture in the 126-strain haploid panel.
Strains are haploid, so there is no dominance. Non-additive variance is epistasis. Two analyses per trait:
 (a) variance components by REML: y = mu + additive (K_A) + additive-by-additive (K_A o K_A) + error; LRT for the epistatic component (boundary-corrected).
 (b) pairwise SNP x SNP interaction scan among the LD-pruned SNPs (r2 < 0.5), LMM-adjusted for relatedness (K_A fitted under the additive-only null)."""
import sys, time
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats, optimize
sys.path.insert(0, "analysis/gwas_dh4148/scripts")
from geno_lib import *
G = R / "gwas"; T = Path("analysis/gwas_dh4148/report/tables")
log = lambda m: print(m, flush=True)
ids, X, sid, names = load_geno(); n = X.shape[1]
P = pd.read_csv(G / "pheno_rint.csv", index_col=0).loc[sid]; traits = list(P.columns)
Xm = np.where(np.isnan(X), np.nanmean(X, axis=1, keepdims=True), X).astype(np.float64)
Z = (Xm - Xm.mean(1, keepdims=True)) / np.maximum(Xm.std(1, keepdims=True), 1e-9)
KA = Z.T @ Z / Z.shape[0]; KA /= np.trace(KA) / n
KAA = KA * KA; KAA /= np.trace(KAA) / n; KAA = KAA - KAA.mean(0, keepdims=True) - KAA.mean(1, keepdims=True) + KAA.mean()   # centred
KAA /= np.trace(KAA) / n
def reml(y, Ks):
    ok = ~np.isnan(y); y = y[ok]; m = ok.sum(); Ks = [K[np.ix_(ok, ok)] for K in Ks]; one = np.ones(m)
    def nll(lp):
        s = np.exp(lp); V = s[-1] * np.eye(m) + sum(si * K for si, K in zip(s[:-1], Ks))
        L = np.linalg.cholesky(V); Vi1 = np.linalg.solve(V, one); Viy = np.linalg.solve(V, y); d = one @ Vi1
        r = y - one * (one @ Viy) / d; return 0.5 * (2 * np.log(np.diag(L)).sum() + np.log(d) + r @ np.linalg.solve(V, r))
    best = None; v = y.var()
    for st in ([0.3, 0.7], [0.1, 0.9], [0.6, 0.4]):
        x0 = np.log(np.array(st[:1] * len(Ks) + [st[1]]) * v / (1 if len(Ks) == 1 else 1) + 1e-6) if len(Ks) == 1 else np.log(np.array([0.3 * v, 0.1 * v, 0.6 * v]))
        r = optimize.minimize(nll, x0, method="Nelder-Mead", options=dict(maxiter=4000, xatol=1e-4, fatol=1e-6))
        if best is None or r.fun < best.fun: best = r
    return np.exp(best.x), -best.fun
rows = []
for t in traits:
    y = P[t].values
    sA, ll1 = reml(y, [KA]); sAE, ll2 = reml(y, [KA, KAA]); lr = max(0.0, 2 * (ll2 - ll1)); p = 0.5 * stats.chi2.sf(lr, 1) if lr > 0 else 1.0
    tot1 = sA.sum(); tot2 = sAE.sum()
    rows.append(dict(trait=t, n=int((~np.isnan(y)).sum()), h2_additive_only=sA[0] / tot1, h2_additive=sAE[0] / tot2, h2_epistatic_AA=sAE[1] / tot2, residual_share=sAE[2] / tot2, LRT=lr, p_epistatic=p))
vc = pd.DataFrame(rows); vc.to_csv(T / "gwas_variance_components_additive_epistatic.csv", index=False); log(vc.round(3).to_string(index=False))
# ---- (b) pairwise scan on pruned SNPs
pr = [l.strip() for l in open(GENO / "pruned_r2_0.5.snps.txt")]; ix = {k: i for i, k in enumerate(ids)}; sel = np.array([ix[k] for k in pr])
H = np.where(np.isnan(X[sel]), np.nan, X[sel]); maj = (np.nanmean(H, axis=1) >= 0.5).astype(float); H = np.where(np.isnan(H), maj[:, None], H)   # missing -> major allele, hard 0/1 calls
p_ = len(sel); log(f"epistasis scan: {p_} pruned SNPs, {p_ * (p_ - 1) // 2:,} pairs")
w_, U = np.linalg.eigh(KA)
Y = P.values.copy(); ok = ~np.isnan(Y)
# null additive fit per trait -> weights in the rotated space (missing phenotypes imputed by the trait mean for the rotation; their weight is not used because traits with missing data use the complete-case subset below)
res = {}
for ti, t in enumerate(traits):
    keep = ok[:, ti]; Kk = KA[np.ix_(keep, keep)]; s, _ = reml(Y[keep, ti], [KA[np.ix_(keep, keep)]]); lam_, Uk = np.linalg.eigh(Kk)
    res[t] = dict(keep=keep, U=Uk, w=1.0 / (s[0] * np.maximum(lam_, 0) + s[1]), y=Uk.T @ Y[keep, ti])
top = {t: [] for t in traits}; allp = {t: [] for t in traits}; t0 = time.time()
Hn = H.T  # n x p
for i in range(p_):
    gi = Hn[:, i]
    for t in traits:
        r = res[t]; k = r["keep"]; Uk, w, y = r["U"], r["w"], r["y"]; sw = np.sqrt(w)
        Gw = sw[:, None] * (Uk.T @ Hn[k, i + 1:]); Iw = sw[:, None] * (Uk.T @ (gi[k][:, None] * Hn[k, i + 1:]))
        if Gw.shape[1] == 0: continue
        B = np.column_stack([sw * (Uk.T @ np.ones(k.sum())), sw * (Uk.T @ gi[k])]); Q, _ = np.linalg.qr(B)
        yr = sw * y; yr = yr - Q @ (Q.T @ yr); Gr = Gw - Q @ (Q.T @ Gw); Ir = Iw - Q @ (Q.T @ Iw)
        a = (Gr * Gr).sum(0); b = (Ir * Gr).sum(0); c = (Ir * Ir).sum(0); tg = Gr.T @ yr; ti_ = Ir.T @ yr; det = a * c - b * b
        # cell counts of the 2x2 genotype classes
        gj = Hn[k, i + 1:]; gk = gi[k][:, None]; n11 = (gk * gj).sum(0); n10 = (gk * (1 - gj)).sum(0); n01 = ((1 - gk) * gj).sum(0); n00 = ((1 - gk) * (1 - gj)).sum(0)
        valid = (det > 1e-8) & (np.minimum.reduce([n11, n10, n01, n00]) >= 5)
        bi = np.where(valid, (a * ti_ - b * tg) / np.where(valid, det, 1), np.nan); bg = np.where(valid, (c * tg - b * ti_) / np.where(valid, det, 1), np.nan)
        rss = yr @ yr - bg * tg - bi * ti_; df = k.sum() - 4; s2 = rss / df; se = np.sqrt(s2 * a / np.where(valid, det, 1)); tt = bi / se; pp = 2 * stats.t.sf(np.abs(tt), df)
        pp = np.where(valid, pp, np.nan); allp[t].append(pp[valid])
        o = np.argsort(np.where(valid, pp, 2))[:5]
        for j in o:
            if valid[j]: top[t].append((pp[j], i, i + 1 + j, bi[j]))
        top[t] = sorted(top[t])[:200]
    if i % 200 == 0: log(f"  i = {i}/{p_}  {time.time() - t0:.0f}s")
out = []; sumr = []
for t in traits:
    ap = np.concatenate(allp[t]); npairs = len(ap); lam = np.median(stats.chi2.isf(ap, 1)) / 0.4549
    sumr.append(dict(trait=t, n_pairs_tested=npairs, bonferroni_threshold=0.05 / npairs, min_p=ap.min(), lambda_gc=lam, n_bonferroni=int((ap < 0.05 / npairs).sum()), n_expected_at_1e_5=npairs * 1e-5, n_at_1e_5=int((ap < 1e-5).sum())))
    for pv, i, j, b in top[t][:20]: out.append(dict(trait=t, snp1=pr[i], snp2=pr[j], p_interaction=pv, interaction_beta_whitened=b))
pd.DataFrame(sumr).to_csv(T / "gwas_epistasis_summary.csv", index=False); pd.DataFrame(out).to_csv(T / "gwas_epistasis_top_pairs.csv", index=False)
log(pd.DataFrame(sumr).round(4).to_string(index=False))
