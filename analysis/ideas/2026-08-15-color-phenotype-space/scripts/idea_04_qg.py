#!/usr/bin/env python3
"""Idea 04 (Quantitative Geneticist): how heritable/repeatable is the screen?

(4a) ICC audit: variance decomposition (strain | plate-well | colony) for a wide
     feature set on Cu=0 late-window colonies; which traits are strain-level
     repeatable; repeatability of the unused within-colony heterogeneity stats.
(4b) Strain x Cu reaction norms for chromatic traits; growth x pigment
     decoupling; Cu x strain interaction F-test.
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
import statsmodels.formula.api as smf
from scipy import stats

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import common as CM

VCFEATS = {
    "Shape_Area": None, "Shape_Compactness": None, "Shape_Circularity": None,
    "Shape_Solidity": None, "Shape_Extent": None,
    "Intensity_MeanIntensity": None, "Intensity_CoefficientVarianceIntensity": None,
    "ColorLab_L*Median": None, "ColorLab_a*Median": None, "ColorLab_b*Median": None,
    "ColorLab_ChromaEstimatedMedian": None, "ColorLab_a*CoeffVar": None,
    "ColorLab_b*CoeffVar": None, "ColorLab_L*CoeffVar": None,
    "ColorHSV_SaturationMedian": None, "ColorHSV_SaturationCoeffVar": None,
    "Colorxy_xMedian": None, "Colorxy_yMedian": None,
    "TextureGray_Entropy-avg-scale05": None,
    "TextureGray_Contrast-avg-scale05": None,
}
HET = ["ColorLab_a*CoeffVar", "ColorLab_b*CoeffVar", "ColorLab_L*CoeffVar",
       "ColorHSV_SaturationCoeffVar", "TextureGray_Entropy-avg-scale05"]


def variance_components(df: pd.DataFrame, feat: str):
    """
    Nested ANOVA (strain -> plate-well -> colony residual).
    Returns variance partition sigma_s (strain), sigma_w (well within strain), sigma_e (colony).
    """
    df = df[["strain_id", "well_id", feat]].dropna().copy()
    df = df[df.strain_id.notna()]
    if len(df) < 50:
        return pd.NA, pd.NA, pd.NA
    df["nv"] = df[feat].astype(float)
    # nice method-of-moments for two-way nested (balanced-ish approximate)
    grand = df.nv.mean()
    S, W = 0.0, 0.0
    s2 = 0.0
    for sid, s in df.groupby("strain_id"):
        S += len(s) * (s.nv.mean() - grand) ** 2
        for wid, w in s.groupby("well_id"):
            W += len(w) * (w.nv.mean() - s.nv.mean()) ** 2
            s2 += ((w.nv - w.nv.mean()) ** 2).sum()
    ms_s, ms_w, ms_e = (S / (df.strain_id.nunique() - 1),
                        W / (df.well_id.nunique() - df.strain_id.nunique()),
                        s2 / (len(df) - df.well_id.nunique()))
    return ms_s, ms_w, ms_e


def icc_components(df, feat):
    ms_s, ms_w, ms_e = variance_components(df, feat)
    if pd.isna(ms_s):
        return np.nan, np.nan, np.nan
    n_wells = df[["strain_id", "well_id"]].drop_duplicates().groupby("strain_id").size()
    n0 = ((len(df) - (df.well_id.nunique() ** 2) / len(df)) /
          (df.strain_id.nunique() - 1)) if df.strain_id.nunique() > 1 else 1
    n1 = ((df.well_id.nunique() - (df.strain_id.nunique() ** 2) / df.strain_id.nunique() * (1 / df.strain_id.nunique())) /
          (df.strain_id.nunique() - 1))
    # simpler: use average cluster size approximation
    n0 = len(df) / df.strain_id.nunique() / 1
    sig_s = ms_s - ms_w
    sig_w = ms_w - ms_e
    sig_e = ms_e
    total = max(sig_s + sig_w + sig_e, 1e-9)
    return max(sig_s, 0) / total, max(sig_w, 0) / total, max(sig_e, 0) / total


def main() -> None:
    print("[idea04] loading extract ...")
    df = CM.read_extract()
    d = df[df.strain_id.notna() & ~df.is_control.astype(bool)
           & (df.tp_h >= 85) & (df.tp_h <= 110)].copy()
    d["well_id"] = (d.run_number.astype(str) + "_" + d.plate_number.astype(str)
                    + "_" + d.well_position.astype(str))
    d["ColorLab_ChromaEstimatedMedian"] = d["ColorLab_ChromaEstimatedMedian"].astype(float)
    print(f"  [idea04] late-window Cu=0 colonies={len(d)}")

    # ---------- 4a ICC audit at Cu=0 ----------
    print("[idea04] ICC audit (Cu=0, strain | well | colony) ...")
    rows = []
    for feat in VCFEATS:
        ms_s, ms_w, ms_e = variance_components(d, feat)
        if pd.isna(ms_s):
            continue
        sig_s, sig_w, sig_e = icc_components(d, feat)
        rows.append({"feature": feat, "ms_strain": float(ms_s), "ms_well": float(ms_w),
                     "ms_colony": float(ms_e), "icc_strain": sig_s,
                     "icc_well": sig_w, "icc_colony": sig_e})
    icc = pd.DataFrame(rows).sort_values("icc_strain", ascending=False)
    CM.save(icc, "idea04_icc_audit.csv")
    het_rows = icc[icc.feature.isin(HET)].copy()
    CM.save(het_rows, "idea04_heterogeneity_repeatability.csv")

    # ---------- 4b reaction norms + GxE ----------
    print("[idea04] reaction norms (well-level, all Cu) ...")
    dn = df[df.strain_id.notna() & ~df.is_control.astype(bool)
            & (df.tp_h >= 85) & (df.tp_h <= 110)].copy()
    dn["well_id"] = (dn.run_number.astype(str) + "_" + dn.plate_number.astype(str)
                     + "_" + dn.well_position.astype(str))
    dn["chr"] = dn["ColorLab_ChromaEstimatedMedian"].astype(float)
    dn["loga"] = np.log(dn["ColorLab_a*Median"].astype(float) + 2.0)
    well = (dn.groupby(["strain_id", "well_id", "copper_mm"])
               .agg(chroma=("chr", "median"),
                    a=("ColorLab_a*Median", "median"),
                    area=("Shape_Area", "median"),
                    species=("species", "first"))
               .reset_index())
    nm = well[well.chroma.notna()].copy()
    rows = []
    for (sid, sp), g in nm.groupby(["strain_id", "species"]):
        g = g.sort_values("copper_mm")
        if g.copper_mm.nunique() < 5:
            continue
        x, y = np.asarray(g.copper_mm), np.log(np.asarray(g.chroma) + 0.5)
        if np.std(x) == 0 or np.std(y) == 0:
            continue
        lin = stats.linregress(x, y)
        rows.append({"strain_id": sid, "species": sp, "n_cu": len(g),
                     "slope_logchroma_per_mM": float(lin.slope),
                     "intercept_logchroma": float(lin.intercept),
                     "r2": float(lin.rvalue ** 2)})
    rn = pd.DataFrame(rows)
    rn = rn.merge(CM.read_meta()[["strain_id", "environment", "origin"]], on="strain_id", how="left")
    CM.save(rn, "idea04_reaction_norms.csv")
    if len(rn) > 10:
        corr = stats.spearmanr(rn.intercept_logchroma, rn.slope_logchroma_per_mM)
        print(f"  [idea04] spearman(baseline chroma, Cu slope) = {corr.statistic:.3f} p={corr.pvalue:.3g}")

    # Cu x strain interaction F (within most-sampled species, well-level)
    print("[idea04] Cu x strain interaction F-test (chroma, R. mucilaginosa) ...")
    muc = nm[nm.species == "Rhodotorula mucilaginosa"].copy()
    muc["logchr"] = np.log(muc.chroma + 0.5)
    if len(muc) > 200:
        full = smf.ols("logchr ~ C(copper_mm) * C(strain_id)", data=muc).fit()
        red = smf.ols("logchr ~ C(copper_mm) + C(strain_id)", data=muc).fit()
        df_int = full.df_resid
        rss_r, rss_f = red.ssr, full.ssr
        f = ((rss_r - rss_f) / (full.df_model - red.df_model)) / (rss_f / df_int)
        p = 1 - stats.f.cdf(f, full.df_model - red.df_model, df_int)
        fx = pd.DataFrame({
            "species": "Rhodotorula mucilaginosa",
            "n_well_rows": len(muc), "n_strains": muc.strain_id.nunique(),
            "F_interaction": f, "p_interaction": p,
            "dof_interaction": int(full.df_model - red.df_model)}, index=[0])
        CM.save(fx, "idea04_gxe_interaction_mucilaginosa.csv")
        print(f"    F={fx.F_interaction.iloc[0]:.1f} p={fx.p_interaction.iloc[0]:.3g}")

    # ---------- figures ----------
    plt.figure(figsize=(8, 5))
    plt.barh(np.arange(len(icc)), icc.icc_strain, color="#1f77b4")
    plt.yticks(np.arange(len(icc)), [f.split("(")[0] if "(" not in f else f for f in icc.feature],
               fontsize=6.5)
    plt.xlabel("ICC_strain (repeatability; strain variance / total)")
    plt.title("Heritability-style audit: strain-level repeatability (Cu=0, late)")
    plt.tight_layout(); plt.savefig(CM.FIGURES / "fig04_icc_audit.png", dpi=150)
    plt.close()

    plt.figure(figsize=(7, 5))
    sp = ["Rhodotorula mucilaginosa", "Rhodotorula paludigena",
          "Rhodotorula sphaerocarpa", "Rhodotorula kratochvilovae"]
    for s in sp:
        d = rn[rn.species == s].dropna()
        if len(d) < 4:
            continue
        x = d.intercept_logchroma
        y = d.slope_logchroma_per_mM
        plt.scatter(x, y, s=10, alpha=0.5, label=s.split()[-1])
    plt.xlabel("baseline log(chroma) @ Cu=0"); plt.ylabel("strain Cu slope of log(chroma)")
    plt.title("Reaction norms: baseline pigment vs Cu sensitivity (strains)")
    plt.legend(fontsize=8); plt.tight_layout()
    plt.savefig(CM.FIGURES / "fig04_reaction_norms.png", dpi=150)
    plt.close()
    print("[idea04] done.")


if __name__ == "__main__":
    sys.exit(main())
