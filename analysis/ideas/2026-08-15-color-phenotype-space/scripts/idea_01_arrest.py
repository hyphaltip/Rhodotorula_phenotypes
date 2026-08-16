#!/usr/bin/env python3
"""Idea 01 (Statistical Physicist): universality of the arrest.

(1a) Master-curve collapse of per-well Shape_Area time-courses across Cu and
     strains; arrest-shape exponent (Weibull shape of the approach to
     saturation) across the Cu field.
(1b) Gibrat / quenched-disorder dispersion: does sd(log Asat) across replicate
     wells widen with Cu?
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

RNG = np.random.default_rng(42)


def well_curves(df: pd.DataFrame) -> pd.DataFrame:
    df = df[df.strain_id.notna() & (df.hours_since_plate_start <= 115)].copy()
    df["well_id"] = (df.run_number.astype(str) + "_" + df.plate_number.astype(str)
                     + "_" + df.well_position.astype(str))
    g = (df.groupby(["well_id", "tp_h"], as_index=False)
           .agg(area=("Shape_Area", "median"),
                copper_mm=("copper_mm", "first"),
                strain_id=("strain_id", "first"),
                species=("species", "first"),
                strain_code=("strain_code", "first")))
    print(f"  [idea01] wells={g.well_id.nunique()} curves rows={len(g)}")
    return g


def collapse_stats(g: pd.DataFrame, cu_keep=(0, 5, 15, 30)):
    rows = []
    for cu in cu_keep:
        sub = g[g.copper_mm == cu]
        res = []
        for wid, w in sub.groupby("well_id"):
            w = w.sort_values("tp_h")
            wp = w[w.area > 0]
            if len(wp) < 8:
                continue
            asat = float(np.quantile(wp.area.values, 0.95))
            if asat <= 0:
                continue
            t = wp.tp_h.values.astype(float)
            a = wp.area.values / asat
            t50 = _t50(t, a)
            if t50 is None:
                res.append({"well_id": wid, "strain_id": w.strain_id.iloc[0],
                            "t50": np.nan, "tau": np.nan, "weibull_k": np.nan,
                            "asat": asat, "reached_half": 0, "n_pass": len(wp)})
                continue
            slope = _slope(t, a, t50)
            tau = 1.0 / slope if slope > 0 else np.nan
            k = _weibull_k(t, a, t50)
            res.append({"well_id": wid, "strain_id": w.strain_id.iloc[0],
                        "t50": t50, "tau": tau, "weibull_k": k,
                        "asat": asat, "reached_half": 1, "n_pass": len(wp)})
        d = pd.DataFrame(res)
        if len(d) == 0:
            continue
        r50 = d.reached_half.mean()
        k_med = d.weibull_k[np.isfinite(d.weibull_k)].to_numpy(float)
        k_samps = [np.nanmedian(k_med[RNG.integers(0, len(k_med), len(k_med))])
                   for _ in range(999)] if len(k_med) else [np.nan]
        rows.append({"copper_mm": cu, "n_wells": len(d),
                     "pct_reached_half": 100 * r50,
                     "median_t50": np.nanmedian(d.t50),
                     "median_tau_h": np.nanmedian(d.tau),
                     "weibull_k_median": float(np.nanmedian(k_med)) if len(k_med) else np.nan,
                     "weibull_k_ci_lo": float(np.percentile(k_samps, 2.5)) if k_samps else np.nan,
                     "weibull_k_ci_hi": float(np.percentile(k_samps, 97.5)) if k_samps else np.nan,
                     })
    return pd.DataFrame(rows)


def _t50(t, a):
    m = np.isfinite(a) & (a > 0)
    if m.sum() < 2:
        return None
    t, a = t[m], a[m]
    if a.max() < 0.5 or a.min() > 0.5:
        return None
    idx = np.where(np.diff(np.sign(a - 0.5)) != 0)[0]
    for i in idx:
        t0, t1 = t[i], t[i + 1]
        a0, a1 = a[i], a[i + 1]
        if a1 == a0:
            continue
        return float(t0 + (0.5 - a0) * (t1 - t0) / (a1 - a0))
    return None


def _slope(t, a, t50):
    # local slope at t50 (central difference over neighbours)
    j = int(np.argmin(np.abs(t - t50)))
    i = max(j - 1, 0)
    k = min(j + 1, len(t) - 1)
    if k == i or t[k] == t[i]:
        return np.nan
    return (a[k] - a[i]) / (t[k] - t[i])


def _weibull_k(t, a, t50):
    # log(-log(1-x)) = k*log(t-t50) + const  for x in (0.5, 0.985)
    tt = np.asarray(t, float)
    aa = np.asarray(a, float)
    m = (tt > t50) & np.isfinite(aa) & (aa > 0.5) & (aa < 0.985)
    if m.sum() < 4:
        return np.nan
    x = np.log(tt[m] - t50)
    y = np.log(-np.log(1.0 - aa[m]))
    res = stats.linregress(x, y)
    return float(res.slope)


def gibrat(g: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (sid, cu), w in g.groupby(["strain_id", "copper_mm"]):
        a = w[w.area > 0].groupby("well_id").area.quantile(0.95)
        a = a[a > 0]
        if len(a) < 2:
            continue
        loga = np.log10(a.values)
        rows.append({"strain_id": sid, "species": w.species.dropna().iloc[0] if w.species.notna().any() else np.nan,
                     "copper_mm": cu, "n_wells": len(a),
                     "sd_log10_asat": float(np.std(loga, ddof=1)),
                     "median_asat": float(np.median(a))})
    return pd.DataFrame(rows)


def main() -> None:
    print("[idea01] loading extract ...")
    df = C.read_extract()
    g = well_curves(df)
    print("[idea01] collapse + arrest exponent ...")
    col = collapse_stats(g)
    C.save(col, "idea01_arrest_collapse.csv")
    print("[idea01] gibrat dispersion ...")
    dis = gibrat(g)
    C.save(dis, "idea01_dispersion.csv")

    # ---- figures ----
    # master curve collapse (mean+SE of x vs (t-t50)/tau, well-level medians)
    plt.figure(figsize=(7, 5))
    for cu in (0.0, 5.0, 15.0, 30.0):
        sub = g[g.copper_mm == cu]
        xs, us = [], []
        for wid, w in sub.groupby("well_id"):
            w = w.sort_values("tp_h")
            wp = w[w.area > 0]
            if len(wp) < 8:
                continue
            asat = float(np.quantile(wp.area.values, 0.95))
            if asat <= 0:
                continue
            nn = np.isfinite(wp.area)
            t, a = wp.tp_h.values[nn].astype(float), wp.area.values[nn] / asat
            t50 = _t50(t, a)
            if t50 is None:
                continue
            u = (t - t50) / max(_slope(t, a, t50), 1e-9)
            m = (u > -1) & (u < 6)
            xs.extend(a[m]); us.extend(u[m])
        us = np.array(us); xs = np.array(xs)
        order = np.argsort(us)
        keep = us[order]; xv = xs[order]
        b = np.linspace(-1, 6, 29)
        mid = (b[:-1] + b[1:]) / 2
        mean, se = [], []
        for i in range(len(b) - 1):
            w = xv[(keep >= b[i]) & (keep < b[i + 1])]
            mean.append(w.mean() if len(w) else np.nan)
            se.append(w.std(ddof=1) / np.sqrt(len(w)) if len(w) else np.nan)
        ok = np.isfinite(mean)
        plt.plot(mid[ok], np.array(mean)[ok], label=f"Cu={cu:g} mM")
        plt.fill_between(mid[ok], np.array(mean)[ok] - 2 * np.array(se)[ok],
                         np.array(mean)[ok] + 2 * np.array(se)[ok], alpha=0.25)
    plt.axhline(0.5, color="grey", ls="--", lw=0.8)
    plt.xlabel("rescaled time  (t-t50)/tau"); plt.ylabel("normalized area  A/A_sat")
    plt.title("Master-curve collapse attempt per Cu (well-level, mean±2SE)")
    plt.legend(); plt.tight_layout()
    plt.savefig(C.FIGURES / "fig01_master_collapse.png", dpi=150)
    plt.close()

    # arrest exponent per Cu boxplot
    plt.figure(figsize=(6.5, 4.5))
    dk = col[["copper_mm", "weibull_k_median", "weibull_k_ci_lo", "weibull_k_ci_hi"]]
    x = np.arange(len(dk))
    plt.errorbar(x, dk.weibull_k_median,
                 yerr=[dk.weibull_k_median - dk.weibull_k_ci_lo,
                       dk.weibull_k_ci_hi - dk.weibull_k_median],
                 fmt="o", capsize=4, ms=6)
    plt.xticks(x, [f"{v:g}" for v in dk.copper_mm])
    plt.xlabel("copper (mM)"); plt.ylabel("Weibull shape k (arrest exponent)")
    plt.title("Arrest-shape exponent vs Cu (bootstrap 95% CI)")
    plt.tight_layout(); plt.savefig(C.FIGURES / "fig01_arrest_exponent.png", dpi=150)
    plt.close()

    # gibrat dispersion vs Cu
    plt.figure(figsize=(7, 5))
    top = dis.species.value_counts().index[:5]
    for sp in top:
        d = dis[dis.species == sp].dropna()
        med = d.sd_log10_asat.groupby(d.copper_mm).median().sort_index()
        plt.plot(med.index, med.values, marker="o", ms=4, label=sp.split()[-1])
    d = dis.dropna()
    if len(d) > 5:
        lin = stats.linregress(d.copper_mm, d.sd_log10_asat)
        xx = np.linspace(d.copper_mm.min(), d.copper_mm.max(), 50)
        plt.plot(xx, lin.intercept + lin.slope * xx, "k--", lw=1.2,
                 label=f"all strains: slope={lin.slope:+.4f}, p={lin.pvalue:.3g}")
    plt.xlabel("copper (mM)"); plt.ylabel("sd(log10 Asat) across replicate wells")
    plt.title("Quenched-disorder dispersion vs Cu (Gibrat axis)")
    plt.legend(); plt.tight_layout(); plt.savefig(C.FIGURES / "fig01_gibrat_dispersion.png", dpi=150)
    plt.close()

    print("[idea01] figures written.")
    print(C.save and "\n[idea01] done.")


if __name__ == "__main__":
    sys.exit(main())
