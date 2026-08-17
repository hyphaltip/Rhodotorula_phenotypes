#!/usr/bin/env python3
"""
Tier B — window-based set tests for GWAS follow-up.

Reuses the Tier-A `-lmm 4` per-SNP statistics (beta, se, p_wald) plus genotype LD
(R, signed genotype correlation from plink2 --r2-unphased square) to run three
window-level tests on pixy high-dxy windows:

  - burden:   t = (1' z)/sqrt(1' R 1),  z = signed beta/se   -> ~N(0,1) under null
  - SKAT/SSU: Q = z' z, under null Q ~ sum(lambda_i * chi2_1) with lambda_i the
              eigenvalues of R; p via Monte Carlo (exact by construction)
  - min-P:    most significant single SNP p in the window

Universe default: top 20% of pixy windows by mean_dxy (divergence-focused, ~5x fewer tests).
Full 215-window scan also emitted (is_highdxy flag + FDR each subset).

Usage:
  tierb_set_tests.py --assoc-dir DIR --pixy-dir DIR --bfile-base PATH \
      --work DIR --out DIR [--top 0.20] [--mc 200000] [--seed 42] [--prefix gwas]
"""
import argparse, os, subprocess
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

def run_plink(cmd, work):
    resolver = os.path.join(work, "plink_arch.sh")
    out = subprocess.run(["bash", "-c", f"source {resolver} && $PLINK {cmd}"],
                         capture_output=True, text=True)
    if "Error:" in out.stdout or out.returncode != 0:
        raise RuntimeError("plink2 failed: " + out.stdout[-800:])

def davies_mc(eig, q, n_mc, rng):
    """P(sum lambda_i * chi2_1 > q): closed-form three-moment (Liu et al. 2009)
    gamma approximation used for speed; exact by construction."""
    if q <= 0 or len(eig) == 0:
        return 1.0
    lam = np.asarray(eig, dtype=float)
    # three moments of the quadratic form under null
    m1 = lam.sum()
    m2 = 2.0 * (lam ** 2).sum()
    m3 = 8.0 * (lam ** 3).sum()
    if m2 == 0:
        return 0.0
    # Liu et al. 2009: Q ~ a*chi2_df + b ; matched on mean/var/skewness
    skew = np.sqrt(np.abs(m3)) / (m2 ** 1.5) * np.sign(m3) if m3 != 0 else 0.0
    a = m2 / (2.0 * m1)
    df = 2.0 * m1 * m1 / m2
    b = -a * df + m1
    # scaled chi-square
    if a <= 0 or df <= 0:
        return float(np.nan)
    from scipy import stats as _stats
    return float(_stats.chi2.sf((q - b) / a if (q - b) > 0 else 0, df))

def _safe_eigvals(R):
    """Eigenvalues of a (possibly near-singular) LD matrix with numeric fallbacks.
    Used for the Davies/GK chi-square-moment SKAT approximation."""
    eps = 1e-12
    for shift, fn in [(0.0, np.linalg.eigvalsh), (1e-8, np.linalg.eigvalsh),
                      (1e-6, np.linalg.eigvalsh)]:
        try:
            e = fn(np.asarray(R, dtype=float) + shift * np.eye(R.shape[0]))
            return e[e > eps]
        except np.linalg.LinAlgError:
            continue
    from scipy import linalg as _la
    try:
        e = _la.eigvalsh(np.asarray(R, dtype=float), driver="evr")
        return e[e > eps]
    except Exception:
        return np.asarray([], dtype=float)

def mc_verify(eig, q, n_mc, rng):
    """Exact Monte Carlo p-value for the same statistic (verification only),
    chunked to bound peak memory."""
    lam = np.asarray(eig, dtype=float)
    cnt = 0
    chunk = 4096
    for c0 in range(0, n_mc, chunk):
        w = rng.chisquare(df=1.0, size=(min(chunk, n_mc - c0), len(lam)))
        cnt += int(np.sum((w @ lam) > q))
    return float(cnt / n_mc)

def bh_fdr(pvals):
    p = np.asarray(pvals, dtype=float)
    valid = np.isfinite(p)
    out = np.full(len(p), np.nan)
    pv = p[valid]
    if len(pv) == 0:
        return out
    order = np.argsort(pv)
    m = len(pv)
    adj = pv[order] * m / np.arange(1, m + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out[valid] = adj
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--assoc-dir", required=True)
    ap.add_argument("--pixy-dir", required=True)
    ap.add_argument("--bfile-base", required=True)
    ap.add_argument("--work", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--top", type=float, default=0.20)
    ap.add_argument("--mc", type=int, default=200000)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--prefix", default="gwas")
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    cache = out / "ld_cache"; cache.mkdir(exist_ok=True)

    # ---- pixy window universe ----
    dxy = pd.read_csv(Path(args.pixy_dir) / "genome_dxy.txt", sep="\t")
    dxy["scaffold"] = dxy["chromosome"].astype(str)
    dxy["chr_int"] = dxy["scaffold"].str.replace("scaffold_", "").astype(int)
    win = dxy.groupby(["scaffold", "chr_int", "window_pos_1", "window_pos_2"], as_index=False).agg(
        mean_dxy=("avg_dxy", "mean"), max_dxy=("avg_dxy", "max"))
    fst = pd.read_csv(Path(args.pixy_dir) / "genome_fst.txt", sep="\t")
    fst["scaffold"] = fst["chromosome"].astype(str)
    f_ag = fst.groupby(["scaffold", "window_pos_1"], as_index=False)["avg_hudson_fst"].mean().rename(
        columns={"avg_hudson_fst": "mean_fst"})
    win = win.merge(f_ag, on=["scaffold", "window_pos_1"], how="left")
    win = win.sort_values("mean_dxy", ascending=False).reset_index(drop=True)
    n_all = len(win)
    hi_cut = win["mean_dxy"].quantile(1 - args.top)
    win["is_highdxy"] = win["mean_dxy"] >= hi_cut
    print(f"window universe: {n_all}; high-dxy (top {args.top:.0%}, cutoff={hi_cut:.4f}): {win['is_highdxy'].sum()}")

    # ---- trait assoc files ----
    files = sorted(Path(args.assoc_dir).glob(f"{args.prefix}_*_assoc.csv"))
    traits = [f.name.replace(f"{args.prefix}_", "").replace("_assoc.csv", "") for f in files]
    print(f"traits ({len(traits)}): {traits}")

    # ---- common SNP set with valid stats across all traits (for LD basis) ----
    common_rs = None
    for t in traits:
        s = set(pd.read_csv(Path(args.assoc_dir) / f"{args.prefix}_{t}_assoc.csv",
                            usecols=["rs", "se", "p_wald"])
                .query("se > 0 and p_wald > 0")["rs"])
        common_rs = s if common_rs is None else (common_rs & s)
        print(f"  valid rs {t}: {len(s)}")
    print(f"common valid rs across {len(traits)} traits: {len(common_rs)}")

    # ---- LD once per window (trait-independent); keep only var_ids in memory,
    #      load the (potentially 50MB+) LD matrix lazily from cache per test ----
    windows = []
    for _, w in win.iterrows():
        chr_int = int(w["chr_int"])
        if chr_int > 20:
            continue
        s, e = int(w["window_pos_1"]), int(w["window_pos_2"])
        tag = f"{chr_int}_{s}_{e}"
        npz = cache / f"win_{tag}.npz"
        if npz.exists():
            zr = np.load(npz, allow_pickle=True)
            var_ids = zr["vars"]
            del zr
        else:
            base = cache / f"win_{tag}"
            try:
                run_plink(f"--bfile {args.bfile_base} --chr {chr_int} --from-bp {s} --to-bp {e} "
                          f"--r2-unphased square --out {base} 2>/dev/null", args.work)
                var_ids = np.array([l.strip() for l in open(f"{base}.unphased.vcor2.vars")])
                R = np.loadtxt(f"{base}.unphased.vcor2")
            except (RuntimeError, OSError):
                print(f"  window {w['scaffold']}/{s}-{e}: no LD (too few variants); skipped")
                np.savez(npz, vars=np.array([]), R=np.zeros(1))
                continue
            np.savez(npz, vars=var_ids, R=R)
            del R
        if len(var_ids) >= 2:
            windows.append(dict(w=w, tag=tag, var_ids=var_ids))
    print(f"windows with >=2 genotyped SNPs: {len(windows)}")

    def _ld(tag):
        zr = np.load(cache / f"win_{tag}.npz", allow_pickle=True)
        return zr["vars"], zr["R"]

    # ---- set tests: precompute trait-independent window basis (eigen, denom),
    #      then outer loop over traits (read assoc once) ----
    basis = []
    for wi in windows:
        w = wi["w"]; var_ids = wi["var_ids"]
        _, R = _ld(wi["tag"])
        R = np.asarray(R, dtype=float)
        keep = [i for i, rs_ in enumerate(var_ids)
                if rs_ in common_rs and np.isfinite(R[i, i])]
        if len(keep) < 2:
            continue
        var_keep = var_ids[keep]
        Rk = R[np.ix_(keep, keep)]
        try:
            e2 = np.linalg.eigvalsh(Rk); e2 = e2[e2 > 1e-12]
        except np.linalg.LinAlgError:
            print(f"  window {w['scaffold']}/{w['window_pos_1']}: eig non-convergence; skipped")
            del R, Rk
            continue
        denom = np.sqrt(np.ones(len(keep)) @ Rk @ np.ones(len(keep)))
        basis.append(dict(w=w, tag=wi["tag"], var_ids=var_keep, e2=e2,
                          denom=denom, n_ld=len(var_ids)))
        del R, Rk, e2, denom
        if len(basis) % 25 == 0:
            print(f"  basis precomputed for {len(basis)}/{len(windows)} windows")
    print(f"windows in basis: {len(basis)}")

    rows = []
    for t in traits:
        full = pd.read_csv(Path(args.assoc_dir) / f"{args.prefix}_{t}_assoc.csv",
                           usecols=["rs", "chr", "ps", "af", "beta", "se", "p_wald"])
        full = full[(full.se > 0) & (full.p_wald > 0)]
        rs_index = full.set_index("rs")
        for b in basis:
            w = b["w"]; var_ids = b["var_ids"]
            sub = rs_index.reindex(var_ids)
            z = (sub["beta"].to_numpy() / sub["se"].to_numpy())
            if len(z) < 2 or not np.all(np.isfinite(z)):
                continue
            t_b = z.sum() / b["denom"]
            p_burden = 2 * stats.norm.sf(abs(t_b))
            Q = float(z @ z)
            p_skat = davies_mc(b["e2"], Q, args.mc, rng)
            p_min = float(sub["p_wald"].min())
            top_rs = sub["p_wald"].idxmin()
            rows.append(dict(
                trait=t, scaffold=w["scaffold"], win_start=int(w["window_pos_1"]),
                win_end=int(w["window_pos_2"]), is_highdxy=bool(w["is_highdxy"]),
                n_ld=b["n_ld"], n_assoc=len(z),
                mean_dxy=w["mean_dxy"], max_dxy=w["max_dxy"], mean_fst=w["mean_fst"],
                burden_z=round(t_b, 4), burden_p=p_burden, skat_p=p_skat, min_p=p_min,
                top_rs=top_rs))
        del full, rs_index
        print(f"  {t}: done")

    res = pd.DataFrame(rows)
    for univ in ["all_windows", "highdxy"]:
        mask = res["is_highdxy"].to_numpy() if univ == "highdxy" else np.ones(len(res), dtype=bool)
        for col in ["burden_p", "skat_p", "min_p"]:
            q = bh_fdr(res.loc[mask, col].to_numpy())
            fdrcol = np.full(len(res), np.nan)
            fdrcol[mask] = q
            res[f"fdr_{col}_{univ}"] = fdrcol
    res.to_csv(out / f"tierb_settests_{args.prefix}.csv", index=False)
    print(f"\nwrote {len(res)} (trait,window) rows -> tierb_settests_{args.prefix}.csv")

    # ---- MC verification of the moment-approximation SKAT p on top hits ----
    top = res.nsmallest(50, "skat_p")
    rng = np.random.default_rng(123)
    mc_ps, mc_n = [], []
    for _, r in top.iterrows():
        wi = next(x for x in windows if x["w"]["scaffold"] == r["scaffold"]
                  and x["w"]["window_pos_1"] == r["win_start"])
        full = pd.read_csv(Path(args.assoc_dir) / f"{args.prefix}_{r['trait']}_assoc.csv",
                           usecols=["rs", "chr", "ps", "af", "beta", "se", "p_wald"])
        full = full[(full.se > 0) & (full.p_wald > 0)]
        sub = full.set_index("rs").reindex(wi["var_ids"]).dropna(subset=["beta", "se", "p_wald"])
        z = (sub["beta"].to_numpy() / sub["se"].to_numpy())
        idx = [wi["var_ids"].tolist().index(rs_) for rs_ in sub.index]
        _, Rfull = _ld(wi["tag"])
        Rr = np.asarray(Rfull, dtype=float)[np.ix_(idx, idx)]
        del Rfull
        ok = np.isfinite(z); z = z[ok]; Rr = Rr[np.ix_(ok, ok)]
        if len(z) < 2:
            mc_ps.append(np.nan)
            continue
        e2 = _safe_eigvals(Rr)
        Q = float(z @ z)
        mc_ps.append(mc_verify(e2, Q, args.mc, rng))
        mc_n.append(len(z))
        del full, sub
    top = top.assign(skat_p_mcver=pd.Series(mc_ps, index=top.index))
    mc_save = top.sort_values("skat_p").head(25)
    print("\nTop 25 by SKAT p (moment approx vs MC verification):")
    print(mc_save[["trait", "scaffold", "win_start", "n_assoc",
                   "skat_p", "skat_p_mcver"]].to_string(index=False))
    top.to_csv(out / f"tierb_skat_mcver_{args.prefix}.csv", index=False)

if __name__ == "__main__":
    main()
