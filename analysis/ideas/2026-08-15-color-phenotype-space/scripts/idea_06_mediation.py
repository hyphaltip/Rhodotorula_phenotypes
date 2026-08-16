#!/usr/bin/env python3
"""Idea 06 (Causal Endpoint Mediation): is Cu-induced pigment loss just a
consequence of growth arrest, or an independent pathway?

Regression-based mediation decomposition on well-level late-window data:
   total  : log(chroma) ~ Cu                  (TE)
   direct : log(chroma) ~ Cu + log(area)      (DE: effect after controlling growth)
Fold-change-style mediation fraction = 1 - b_direct/b_total. Bootstrap CIs.
Checked for chroma and for a* (redness). Species as covariate.
"""
from __future__ import annotations

import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import common as CM


def fit_models(d: pd.DataFrame, yvar: str, species_fix: bool = True):
    d = d[[yvar, "logical_area", "copper_mm", "species"]].dropna()
    d = d[d.species.notna()]
    if len(d) < 100:
        return None
    d = d.copy()
    d["logarea"] = np.log(np.asarray(d.logical_area, dtype=float) + 1.0)
    d["y"] = np.asarray(d[yvar], dtype=float)
    d["copper_mm"] = np.asarray(d.copper_mm, dtype=float)
    # total: linear in Cu (well-level; single linear slope is the comparable summary)
    tt = sm.OLS(d.y, sm.add_constant(d.copper_mm)).fit()
    bt = float(tt.params["copper_mm"]); pt = float(tt.pvalues["copper_mm"])
    # direct: Cu effect after controlling growth (and species)
    if species_fix and d.species.nunique() > 1:
        ex = pd.concat([d[["copper_mm", "logarea"]],
                        pd.get_dummies(d.species, prefix="sp", drop_first=True)], axis=1)
        ex = ex.astype(float)
        dd = sm.OLS(d.y, sm.add_constant(ex)).fit()
    else:
        dd = sm.OLS(d.y, sm.add_constant(d[["copper_mm", "logarea"]].astype(float))).fit()
    bd, pdv = float(dd.params["copper_mm"]), float(dd.pvalues["copper_mm"])
    return {"n": len(d), "b_total": bt, "p_total": pt,
            "b_direct_growth": bd, "p_direct": pdv,
            "r2_total": float(tt.rsquared), "r2_direct": float(dd.rsquared)}


def main() -> None:
    print("[idea06] loading extract ...")
    df = CM.read_extract()
    d = df[df.strain_id.notna() & ~df.is_control.astype(bool)
           & (df.tp_h >= 85) & (df.tp_h <= 110)].copy()
    d["chr"] = np.log(d["ColorLab_ChromaEstimatedMedian"].astype(float) + 0.5)
    d["a"] = d["ColorLab_a*Median"].astype(float)
    d["area"] = d["Shape_Area"].astype(float)
    well = (d.groupby(["strain_id", "well_position"])
              .agg(chroma=("chr", "median"), a=("a", "median"),
                   area=("area", "median"), copper_mm=("copper_mm", "first"),
                   species=("species", "first")).reset_index())
    well = well.dropna().copy()
    well["logical_area"] = well.area

    results = {}
    for name, yv in [("chroma", "chroma"), ("redness_a", "a")]:
        r = fit_models(well, yv)
        if r:
            r["mediation_fraction"] = 1 - r["b_direct_growth"] / r["b_total"]
            results[name] = r
            print(f"  [idea06] {name}: n={r['n']} b_total={r['b_total']:.4f} "
                  f"b_direct={r['b_direct_growth']:.4f} frac={r['mediation_fraction']:.2f}")

    # bootstrap CIs for mediation fraction (chroma)
    bs = []
    rng = np.random.default_rng(42)
    for _ in range(500):
        w = well.sample(frac=1.0, replace=True, random_state=rng)
        r = fit_models(w, "chroma")
        if r and r["b_total"] != 0:
            bs.append(1 - r["b_direct_growth"] / r["b_total"])
    bs = np.array(bs)
    lo, hi = np.percentile(bs, [2.5, 97.5])
    print(f"  [idea06] boot CI mediation fraction(95%): [{lo:.2f}, {hi:.2f}]")
    results["chroma"]["boot_ci_lo"] = lo
    results["chroma"]["boot_ci_hi"] = hi

    out = pd.DataFrame(results).T.reset_index().rename(columns={"index": "pathway"})
    CM.save(out, "idea06_mediation_decomposition.csv")

    # figure: binned log(chroma) vs log(area) colored by Cu (decoupling check)
    dbin = well.sample(n=min(len(well), 6000), random_state=1)
    plt.figure(figsize=(6.5, 5))
    sc = plt.scatter(np.log(dbin.logical_area + 1), dbin.chroma, c=dbin.copper_mm,
                     s=6, cmap="viridis_r", alpha=0.55)
    plt.colorbar(sc, label="Cu (mM)")
    plt.xlabel("log colony area"); plt.ylabel("log chroma")
    plt.title("Growth vs pigment with Cu (decoupling check)")
    plt.tight_layout(); plt.savefig(CM.FIGURES / "fig06_growth_pigment_decouple.png", dpi=150)
    plt.close()
    print("[idea06] done.")


if __name__ == "__main__":
    sys.exit(main())
