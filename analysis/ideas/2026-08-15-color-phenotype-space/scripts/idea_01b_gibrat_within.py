#!/usr/bin/env python3
"""Idea 01b (Statistical Physicist) — resolve the Gibrat dispersion ambiguity.

The v1 gibrat figure pooled all (strain x Cu) rows into one linear fit of
sd(log10 Asat) ~ Cu. That is ambiguous because (i) R. mucilaginosa is ~2/3 of
the rows, so the all-strain slope may be composition-driven, and (ii) sd(log10)
is strongly anti-correlated with colony size (r2 ~ -0.6): Cu shrinks colonies,
so a rising dispersion may be a size/floor artifact, not a widening of
quenched disorder per genotype.

This script separates the two confounds:
  1. WITHIN-STRAIN repeated-measure slopes: for each strain, does its own
     across-well dispersion widen with Cu? (immune to between-strain composition)
  2. SIZE-CONTROLLED slopes: re-fit sd ~ Cu | log10(median Asat) per strain; the
     residual Cu slope tests whether any widening survives the size coupling.
  3. Species x Cu: do within-strain slopes differ among species?
"""
from __future__ import annotations

import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import common as C

RNG = np.random.default_rng(7)
MIN_CU_LEVELS = 4      # Cu levels a strain must span for a within-strain slope
MIN_WELLS = 3          # wells per strain x Cu used to estimate sd(log10)


def within_strain_slopes(dis: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (sid, sp), g in dis.groupby(["strain_id", "species"]):
        g = g.dropna(subset=["copper_mm", "sd_log10_asat", "l10med"])
        g = g[g.n_wells >= MIN_WELLS]
        if len(g) < MIN_CU_LEVELS:
            continue
        cu = g.copper_mm.values.astype(float)
        sd = g.sd_log10_asat.values.astype(float)
        lm1 = stats.linregress(cu, sd)
        lm2 = _partial_two_predictors(cu, sd, g.l10med.values.astype(float))
        rho, p_rho = stats.spearmanr(cu, sd)
        rows.append({
            "strain_id": sid, "species": sp, "n_cu": len(g),
            "cu_range": float(cu.max() - cu.min()),
            "slope_sd_cu": float(lm1.slope),        # raw within-strain slope
            "r_sd_cu": float(lm1.rvalue),
            "partial_slope_sd_cu": lm2[0],          # sd ~ cu | log10 median area
            "partial_p": lm2[1],
            "spearman_rho": float(rho) if not np.isnan(rho) else np.nan,
            "spearman_p": float(p_rho) if not np.isnan(p_rho) else np.nan,
            "l10med_fixed": float(np.mean(g.l10med)),
        })
    return pd.DataFrame(rows)


def _partial_two_predictors(x, y, z):
    # y ~ b0 + bx*x + bz*z  by least squares; returns (bx_coef, bx_pvalue)
    X = np.column_stack([np.ones_like(x), x, z])
    XtX = X.T @ X
    try:
        beta = np.linalg.solve(XtX, X.T @ y)
    except np.linalg.LinAlgError:
        return (np.nan, np.nan)
    resid = y - X @ beta
    dfe = len(y) - X.shape[1]
    if dfe < 1:
        return (beta[1], np.nan)
    sigma2 = resid @ resid / dfe
    cov = np.linalg.inv(XtX) * sigma2
    se = np.sqrt(np.diag(cov)[1])
    if se == 0 or not np.isfinite(se):
        return (beta[1], np.nan)
    t = beta[1] / se
    p = 2 * stats.t.sf(abs(t), dfe)
    return (beta[1], p)


def species_table(ws: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for sp, g in ws.groupby("species"):
        if len(g) < 5:
            continue
        n = len(g)
        def boot_median(v):
            s = [np.median(v[RNG.integers(0, n, n)]) for _ in range(499)]
            return float(np.percentile(s, 2.5)), float(np.percentile(s, 97.5))
        med_raw, (lo1, hi1) = np.median(g.slope_sd_cu), boot_median(g.slope_sd_cu.to_numpy(float))
        med_partial, (lo2, hi2) = np.median(g.partial_slope_sd_cu), boot_median(g.partial_slope_sd_cu.to_numpy(float))
        # one-sample tests on strain-level slopes
        w_raw = stats.wilcoxon(g.slope_sd_cu.to_numpy(float), alternative="two-sided") if n >= 6 else (np.nan, np.nan)
        w_part = stats.wilcoxon(g.partial_slope_sd_cu.to_numpy(float), alternative="two-sided") if n >= 6 else (np.nan, np.nan)
        rows.append({
            "species": sp, "n_strains": n,
            "median_slope_raw": float(med_raw), "raw_ci_lo": lo1, "raw_ci_hi": hi1,
            "pct_positive_raw": float(100 * (g.slope_sd_cu > 0).mean()),
            "wilcoxon_raw_p": float(w_raw[1]) if n >= 6 and not np.isnan(w_raw[1]) else np.nan,
            "median_slope_controlled": float(med_partial), "partial_ci_lo": lo2, "partial_ci_hi": hi2,
            "pct_positive_controlled": float(100 * (g.partial_slope_sd_cu > 0).mean()),
            "wilcoxon_partial_p": float(w_part[1]) if n >= 6 and not np.isnan(w_part[1]) else np.nan,
        })
    return pd.DataFrame(rows)


def size_confound(dis: pd.DataFrame) -> dict:
    d = dis.dropna(subset=["sd_log10_asat", "l10med"])
    rho, p = stats.spearmanr(d.sd_log10_asat, d.l10med)
    # within-species partial Spearman of sd ~ cu controlling median size (rank residual)
    rows = []
    for sp, g in d.groupby("species"):
        if len(g) < 30:
            continue
        g = g.copy()
        g["rk_cu"] = g.copper_mm.rank()
        g["rk_sd"] = g.sd_log10_asat.rank()
        g["rk_m"] = g.l10med.rank()
        # residualize ranks
        coef, rho_p = _partial_two_predictors(g.rk_cu.to_numpy(float), g.rk_sd.to_numpy(float), g.rk_m.to_numpy(float))
        rows.append({"species": sp, "n": len(g), "partial_rho_sd_cu_ctl_size": coef,
                     "partial_rank_p": rho_p if np.isfinite(rho_p) else np.nan})
    return {"overall_rho_sd_size": float(rho), "overall_p": float(p)}, pd.DataFrame(rows)


def main() -> None:
    print("[idea01b] loading dispersion ...")
    dis = pd.read_csv(C.RESULTS / "idea01_dispersion.csv")
    dis["copper_mm"] = pd.to_numeric(dis.copper_mm, errors="coerce")
    dis["l10med"] = np.log10(dis.median_asat.clip(lower=1.0))
    dis = dis[dis.n_wells >= MIN_WELLS].copy()

    print("[idea01b] within-strain slopes ...")
    ws = within_strain_slopes(dis)
    C.save(ws, "idea01b_within_strain.csv")
    print(f"  strains with >= {MIN_CU_LEVELS} Cu levels: {len(ws)}")

    print("[idea01b] species aggregation ...")
    st = species_table(ws)
    C.save(st, "idea01b_species_table.csv")
    print(st.round(4).to_string(index=False))

    print("[idea01b] size confound ...")
    conf, conf_ps = size_confound(dis)
    C.save(conf_ps, "idea01b_size_confound.csv")
    print(f"  overall spearman(sd, log10 median area) = {conf['overall_rho_sd_size']:.3f} (p={conf['overall_p']:.2e})")

    # overall pooled one-sample tests (immune to composition)
    pool = ws.dropna(subset=["slope_sd_cu", "partial_slope_sd_cu"])
    w_raw_all = stats.wilcoxon(pool.slope_sd_cu.to_numpy(float))
    w_part_all = stats.wilcoxon(pool.partial_slope_sd_cu.to_numpy(float))
    print(f"  ALL strains n={len(pool)}: raw median slope={np.median(pool.slope_sd_cu):+.4f} "
          f"wilcoxon p={w_raw_all.pvalue:.2e} | %positive={100*(pool.slope_sd_cu>0).mean():.1f}%")
    print(f"  ALL strains: size-controlled median slope={np.median(pool.partial_slope_sd_cu):+.4f} "
          f"wilcoxon p={w_part_all.pvalue:.2e} | %positive={100*(pool.partial_slope_sd_cu>0).mean():.1f}%")
    kw = stats.kruskal(*[g.partial_slope_sd_cu.to_numpy(float) for _, g in pool.groupby("species") if len(g) >= 5])
    print(f"  species x Cu (Kruskal on within-strain controlled slopes): H={kw.statistic:.2f}, p={kw.pvalue:.2e}")

    # ---- figures ----
    plot_species_slopes(st)
    plot_size_scatter(dis, ws)
    print("[idea01b] figures written.")
    print("[idea01b] done.")


def plot_species_slopes(st: pd.DataFrame) -> None:
    st = st.sort_values("median_slope_controlled")
    x = np.arange(len(st))
    plt.figure(figsize=(8, 5))
    plt.errorbar(x, st.median_slope_controlled,
                 yerr=[st.median_slope_controlled - st.partial_ci_lo,
                       st.partial_ci_hi - st.median_slope_controlled],
                 fmt="o", capsize=4, ms=7)
    plt.axhline(0, color="grey", ls="--", lw=1)
    plt.xticks(x, [s.split()[-1] for s in st.species], rotation=30, ha="right")
    plt.ylabel("within-strain slope  d[sd(log10 Asat)]/d[Cu]  (size-controlled)")
    plt.title("Within-genotype dispersion widening per Cu (median + bootstrap 95% CI);\n"
              "horizontal = no widening once colony size is controlled")
    plt.tight_layout()
    plt.savefig(C.FIGURES / "fig01b_within_strain_slopes.png", dpi=150)
    plt.close()


def plot_size_scatter(dis: pd.DataFrame, ws: pd.DataFrame) -> None:
    d = dis.dropna(subset=["sd_log10_asat", "l10med"]).copy()
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    sc = ax[0].scatter(d.l10med, d.sd_log10_asat, c=d.copper_mm, s=8, alpha=0.5,
                       cmap="viridis")
    cb = fig.colorbar(sc, ax=ax[0]); cb.set_label("copper (mM)")
    ax[0].set_xlabel("log10 median Asat (colony size)"); ax[0].set_ylabel("sd(log10 Asat)")
    ax[0].set_title("Dispersion vs colony size at each Cu\n(pooled strain x Cu rows)")

    sp_order = [s for s in ws.species.value_counts().index if ws[ws.species == s].shape[0] >= 5][:5]
    for sp in sp_order:
        q = ws[ws.species == sp]
        ax[1].hist(q.partial_slope_sd_cu.dropna().to_numpy(float), bins=24,
                   alpha=0.35, label=sp.split()[-1])
    ax[1].axvline(0, color="grey", ls="--", lw=1)
    ax[1].set_xlabel("within-strain slope  d[sd(log10)]/d[Cu]  (size-controlled)")
    ax[1].set_ylabel("strains")
    ax[1].set_title("Distribution of within-strain dispersion slopes")
    ax[1].legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(C.FIGURES / "fig01b_size_vs_dispersion.png", dpi=150)
    plt.close()


if __name__ == "__main__":
    sys.exit(main())
