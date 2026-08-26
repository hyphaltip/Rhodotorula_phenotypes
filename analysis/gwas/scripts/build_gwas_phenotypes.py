#!/usr/bin/env python3
"""GWAS phenotype derivation for the analysis/gwas/ port (213-strain rebuild).

Ported from analysis/ideas/2026-08-15-color-phenotype-space/scripts/build_gwas_phenotypes.py
(original design per 09-NEXT-GWAS-DESIGN.md S6), re-pointed at the rebuilt
213-strain .fam (Task 7) instead of the prior 201-strain run's scratch path.

Builds strain-level traits for R. mucilaginosa from the per-colony time-course
extract, aligned to the rebuilt GEMMA .fam strain order, so each trait is a
single GEMMA scan away. Outputs a wide per-strain CSV plus a fam-order-aligned CSV.

Produced (per 09-NEXT-GWAS-DESIGN.md S6):
  A. Control color (Cu=0, YPD)  : chroma, saturation, brightness  (clone-mean over plates)
  B. Colony size growth rate    : late-window median area (clone-mean over plates)
  C. Copper-response block      : AUC per dose, dose-slope, IC50_est, AUC_ratio,
                                  AUC(30)/AUC(0) resilience  (multi-trait Fisher/TATES target)

Expected use: `pixi run python analysis/gwas/scripts/build_gwas_phenotypes.py`.
"""
from __future__ import annotations

import pathlib
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import common as CM

# -----------------------------------------------------------------------------
# Config
# -----------------------------------------------------------------------------
SPECIES = "Rhodotorula mucilaginosa"
# late-window for "endpoint" color/size at control (used by idea04 for chroma traits)
WIN_LO, WIN_HI = 85.0, 110.0
MIN_PLATES = 2            # min replicate plates/strain for clone-mean reliable
COPPER_DOSES = sorted([0, 5, 10, 15, 20, 25, 30])
REF_CU = 10.0             # reference dose for AUC_ratio
# genotyped strain order from the rebuilt kinship's .fam (col-2) -- Task 7 output.
FAM = pathlib.Path(__file__).resolve().parents[1] / "results/gwas/grm_conditioning/rebuilt_kinship/gwas.pruned.fam"


def main() -> None:
    print("[gwas_pheno] loading extract ...")
    df = CM.read_extract()
    d = df[(df.species == SPECIES)].copy()
    d = d[d.strain_code.notna()]
    d["area"] = d["Shape_Area"].astype(float)
    d["chroma"] = d["ColorLab_ChromaEstimatedMedian"].astype(float)
    d["sat"] = d["ColorHSV_SaturationMedian"].astype(float)
    d["bright"] = d["ColorHSV_BrightnessMedian"].astype(float)
    d["tp_h"] = d["tp_h"].astype(float)
    d["cu"] = d["copper_mm"].astype(float)
    d["pid"] = (d.run_number.astype(str) + "_" + d.plate_number.astype(str)
                + "_" + d.well_position.astype(str))
    print(f"  [gwas_pheno] Rmuc rows={len(d)}, strains={d.strain_code.nunique()}")

    assert FAM.exists(), f"rebuilt .fam not found: {FAM} -- run Task 7's rebuild_tiers_abc.sh first"
    fam_order = [ln.split()[1] for ln in open(FAM).read().splitlines() if ln.split()]
    assert len(fam_order) == len(set(fam_order)), "duplicate strain IDs in rebuilt .fam"

    # ---------------- A. Control color (clone-mean over plates) ----------------
    print("[gwas_pheno] A. control color (Cu=0, late window, clone-mean over plates) ...")
    a = d[(d.cu == 0) & (d.tp_h >= WIN_LO) & (d.tp_h <= WIN_HI)].copy()
    plate = (a.groupby(["strain_code", "run_number", "plate_number"])
               .agg(chroma=("chroma", "median"), sat=("sat", "median"),
                    bright=("bright", "median"), area=("area", "median"))
               .reset_index())
    plate = plate[plate.chroma.notna()]
    nplate = plate.groupby("strain_code").size()
    plate_mean = (plate.groupby("strain_code")
                       .agg(chroma=("chroma", "mean"), sat=("sat", "mean"),
                            bright=("bright", "mean"), area=("area", "mean"))
                       .reset_index())
    plate_mean["n_plate"] = plate_mean.strain_code.map(nplate)
    plate_mean["clone_mean_area"] = np.log10(plate_mean.area)
    plate_mean["n_rep"] = plate_mean.strain_code.map(nplate)
    color = plate_mean[["strain_code", "chroma", "sat", "bright", "clone_mean_area", "n_plate"]].copy()
    print("   color traits ready, strains=", len(color))

    # ---------------- C. Copper-response block ---------------------------------
    print("[gwas_pheno] C. copper-response block (dose-response per strain) ...")
    rows = []
    for (sc, cu), g in d.groupby(["strain_code", "cu"]):
        auc = []
        late = []
        for pid, w in g.groupby("pid"):
            w = w.sort_values("tp_h")
            t = w["tp_h"].to_numpy(dtype=float)
            ar = w["area"].to_numpy(dtype=float)
            if len(t) < 5:
                continue
            auc.append(np.trapezoid(ar, t))
            m = (t >= WIN_LO) & (t <= WIN_HI)
            late.append(np.median(ar[m]) if m.sum() else np.nan)
        rows.append({"strain_code": sc, "cu": cu, "AUC": float(np.mean(auc)) if auc else np.nan,
                     "late_area": float(np.nanmean(late)) if late else np.nan,
                     "n_well": len(auc)})
    dr = pd.DataFrame(rows)
    aucw = dr.pivot_table(index="strain_code", columns="cu", values="AUC")
    cu_resp = pd.DataFrame(index=aucw.index)
    cu_resp["AUC_0"] = aucw.get(0.0)
    cu_resp["AUC_10"] = aucw.get(10.0)
    cu_resp["AUC_20"] = aucw.get(20.0)
    cu_resp["AUC_30"] = aucw.get(30.0)
    cu_resp["AUC_ratio_10"] = cu_resp.AUC_10 / cu_resp.AUC_0
    cu_resp["resilience_30"] = cu_resp.AUC_30 / cu_resp.AUC_0
    cu_resp = cu_resp.reset_index()
    dresp_rows = []
    for sc, g in dr[dr.cu >= 0].groupby("strain_code"):
        g = g.sort_values("cu")
        y = np.log(g.late_area + 1)
        x = np.asarray(g.cu, float)
        if np.std(x) == 0 or np.std(y) == 0 or len(g) < 4:
            continue
        lin = stats.linregress(x, y)
        dresp_rows.append({"strain_code": sc, "cu_dose_slope": float(lin.slope),
                           "cu_dose_r2": float(lin.rvalue ** 2), "n_cu": len(g)})
    ds = pd.DataFrame(dresp_rows)
    cu_resp = cu_resp.merge(ds, on="strain_code", how="left")
    ic_rows = []
    for sc, g in dr[dr.cu >= 0].groupby("strain_code"):
        g = g.sort_values("cu").reset_index(drop=True)
        if 0.0 not in set(g.cu) or np.isnan(g.late_area.iloc[0]) or g.late_area.max() <= 0:
            continue
        base = g.late_area.iloc[0]
        target = 0.5 * base
        ic50 = np.nan
        prev = base
        for i in range(1, len(g)):
            cur = g.late_area.iloc[i]
            if np.isnan(cur):
                continue
            if prev >= target >= cur:
                x0, x1 = g.cu.iloc[i - 1], g.cu.iloc[i]
                ic50 = x0 + (target - prev) * (x1 - x0) / (cur - prev)
                break
            prev = cur
        ic_rows.append({"strain_code": sc, "IC50_est": ic50})
    ic50 = pd.DataFrame(ic_rows)
    cu_resp = cu_resp.merge(ic50, on="strain_code", how="left")
    print("   copper traits ready, strains=", len(cu_resp))

    # ---------------- Merge all and emit ---------------------------------------
    out = color.merge(cu_resp, on="strain_code", how="outer")
    out.to_csv(CM.RESULTS / "gwas_next_phenotypes.csv", index=False)

    # ---------------- Align to GEMMA .fam order --------------------------------
    aligned = out[out.strain_code.isin(fam_order)].copy()
    idx = aligned.strain_code.map({s: i for i, s in enumerate(fam_order)})
    aligned = aligned.loc[idx.sort_values().index].reset_index(drop=True)
    aligned.to_csv(CM.RESULTS / "gwas_next_phenotypes_fam_order.csv", index=False)
    print(f"  [gwas_pheno] aligned {len(aligned)}/{len(fam_order)} strains to gwas.pruned.fam")
    missing = set(fam_order) - set(aligned.strain_code)
    if missing:
        print("   WARN strains in .fam without extract phenotypes:", sorted(missing)[:20])

    print("[gwas_pheno] done.")


if __name__ == "__main__":
    sys.exit(main())
