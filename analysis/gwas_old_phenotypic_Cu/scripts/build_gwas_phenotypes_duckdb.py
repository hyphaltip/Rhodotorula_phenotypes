#!/usr/bin/env python3
"""Re-derive the GWAS strain-level traits from heavy_metal_measurement (Copper).

Port of build_gwas_phenotypes.py to the new DuckDB table. It does NOT overwrite
analysis/gwas/results/gwas_next_phenotypes*.csv. It writes to --out-dir and
compares each trait to the existing (old-data) table.

Rules (all explicit; see analysis/gwas/DUCKDB_RETOOL_PLAN.md):
- Metal = Copper only. Conc 0..30 (same dose grid as old copper_mm).
- Species = 'Rhodotorula mucilaginosa' per strain_info.species (current DB value).
- strain_code = strain_info.sample_name (TFCN_*), the ID in the reviewed match table.
- Hours: capture_datetime minus the earliest capture_datetime of the same (run_number,
  plate_position) -- the old v_image definition (checked in the pre-metal-replace DB:
  min(imaged_at) OVER (PARTITION BY run_number, plate_number)). tp_h = round(hours).
- Well id (pid) = run_number|plate_position|Grid_RowMajorIdx (one colony/well/capture).
- Window for late traits: 85-110 h on tp_h, as in the old script.
- Color: ColorLab_{L*,a*,b*}GeoMedian replaces the old per-colony *Median. chroma =
  sqrt(a^2+b^2) from GeoMedian (old: ColorLab_ChromaEstimatedMedian). sat/bright use
  ColorHSV_SaturationRobustMean / ColorHSV_ValueRobustMean (old: HSV *Median).
  These three are APPROXIMATIONS, not the same statistic; see the comparison output.
- Texture: Texture_<M>-avg-scale05 replaces TextureGray_<M>-avg-scale05.
- Rows with NULL strain_id (Control-N and NULL) are excluded and counted.
"""
from __future__ import annotations
import argparse, pathlib, sys
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import duckdb_inputs as DI

REPO = DI.REPO
SPECIES = "Rhodotorula mucilaginosa"
WIN_LO, WIN_HI = 85.0, 110.0
TEX = ["AngularSecondMoment", "Contrast", "Correlation", "HaralickVariance",
       "InverseDifferenceMoment", "SumAverage", "SumVariance", "SumEntropy",
       "Entropy", "DiffVariance", "DiffEntropy", "InfoCorrelation1", "InfoCorrelation2"]


def load(con) -> pd.DataFrame:
    texsel = ", ".join(f'"Texture_{m}-avg-scale05" AS tex_{m}' for m in TEX)
    sql = f"""
        SELECT strain_id, run_number, plate_position, Grid_RowMajorIdx, capture_datetime,
               Concentration AS cu, Shape_Area AS area,
               "ColorLab_L*GeoMedian" AS lab_L, "ColorLab_a*GeoMedian" AS lab_a,
               "ColorLab_b*GeoMedian" AS lab_b,
               ColorHSV_SaturationRobustMean AS sat, ColorHSV_ValueRobustMean AS bright,
               {texsel}
        FROM heavy_metal_measurement WHERE Metal = 'Copper'"""
    return con.execute(sql).df()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()
    out = pathlib.Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)
    con = DI.connect_readonly()
    d = load(con)
    si = con.execute("SELECT strain_id, sample_name AS strain_code, species FROM strain_info WHERE NOT is_control").df()
    print(f"[pheno] Copper rows={len(d)}")
    assert len(d) == 155470, f"expected 155,470 Copper rows, got {len(d)}"
    n_null = d.strain_id.isna().sum()
    d = d.merge(si, on="strain_id", how="left")
    n_ctrl_or_null = d.strain_code.isna().sum()
    print(f"  rows with NULL strain_id={n_null}; rows without strain_info match (NULL or Control)={n_ctrl_or_null}")
    d = d[d.strain_code.notna()]
    d = d[d.species == SPECIES].copy()
    print(f"  after species=={SPECIES!r}: rows={len(d)}, strains(strain_code)={d.strain_code.nunique()}, strain_ids={d.strain_id.nunique()}")
    # plate start must come from ALL Copper rows of the plate, not only retained ones
    allt0 = con.execute("SELECT run_number, plate_position, min(capture_datetime) t0 FROM heavy_metal_measurement WHERE Metal='Copper' GROUP BY 1,2").df()
    d = d.merge(allt0, on=["run_number", "plate_position"], how="left")
    d["hours"] = (d.capture_datetime - d.t0).dt.total_seconds() / 3600.0
    assert d.hours.notna().all() and d.hours.min() >= 0
    print(f"  hours range {d.hours.min():.2f}-{d.hours.max():.2f}")
    d["tp_h"] = d.hours.round()
    d["pid"] = d.run_number + "|" + d.plate_position.astype(str) + "|" + d.Grid_RowMajorIdx.astype(str)
    dup = d.duplicated(["pid", "tp_h"]).sum()
    print(f"  duplicate (well, rounded hour) rows={dup} -> collapsed to the median per (well, rounded hour)")
    d = (d.groupby(["strain_code", "strain_id", "cu", "pid", "run_number", "plate_position", "tp_h"], as_index=False)
           [["area", "lab_L", "lab_a", "lab_b", "sat", "bright"] + [f"tex_{m}" for m in TEX]].median())
    print(f"  rows after collapsing: {len(d)}")
    d["chroma"] = np.hypot(d.lab_a, d.lab_b)

    a = d[(d.cu == 0) & (d.tp_h >= WIN_LO) & (d.tp_h <= WIN_HI)].copy()
    print(f"  control (Cu=0) rows in {WIN_LO}-{WIN_HI} h: {len(a)}")
    tex_cols = [f"tex_{m}" for m in TEX]
    cols = ["chroma", "sat", "bright", "area", "lab_L", "lab_a", "lab_b"] + tex_cols
    plate = a.groupby(["strain_code", "run_number", "plate_position"])[cols].median().reset_index()
    nplate = plate.groupby("strain_code").size()
    pm = plate.groupby("strain_code")[cols].mean().reset_index()
    pm["clone_mean_area"] = np.log10(pm.area)
    pm["n_plate"] = pm.strain_code.map(nplate)
    color = pm.drop(columns=tex_cols + ["area"])
    texture = pm[["strain_code"] + tex_cols]

    rows = []
    for (sc, cu), g in d.groupby(["strain_code", "cu"]):
        auc, late = [], []
        for pid, w in g.groupby("pid"):
            w = w.sort_values("tp_h")
            t = w.tp_h.to_numpy(float); ar = w.area.to_numpy(float)
            if len(t) < 5:
                continue
            auc.append(np.trapezoid(ar, t))
            m = (t >= WIN_LO) & (t <= WIN_HI)
            late.append(np.median(ar[m]) if m.sum() else np.nan)
        rows.append({"strain_code": sc, "cu": cu, "AUC": float(np.mean(auc)) if auc else np.nan,
                     "late_area": float(np.nanmean(late)) if late and not np.all(np.isnan(late)) else np.nan,
                     "n_well": len(auc)})
    dr = pd.DataFrame(rows)
    aucw = dr.pivot_table(index="strain_code", columns="cu", values="AUC")
    cr = pd.DataFrame(index=aucw.index)
    for k in (0, 10, 20, 30):
        cr[f"AUC_{k}"] = aucw.get(float(k))
    cr["AUC_ratio_10"] = cr.AUC_10 / cr.AUC_0
    cr["resilience_30"] = cr.AUC_30 / cr.AUC_0
    cr = cr.reset_index()
    sl = []
    for sc, g in dr.groupby("strain_code"):
        g = g.sort_values("cu").dropna(subset=["late_area"])
        y = np.log(g.late_area + 1); x = g.cu.to_numpy(float)
        if len(g) < 4 or np.std(x) == 0 or np.std(y) == 0:
            continue
        lin = stats.linregress(x, y)
        sl.append({"strain_code": sc, "cu_dose_slope": float(lin.slope), "cu_dose_r2": float(lin.rvalue ** 2), "n_cu": len(g)})
    cr = cr.merge(pd.DataFrame(sl), on="strain_code", how="left")
    ic = []
    for sc, g in dr.groupby("strain_code"):
        g = g.sort_values("cu").reset_index(drop=True)
        if 0.0 not in set(g.cu) or np.isnan(g.late_area.iloc[0]) or g.late_area.max() <= 0:
            continue
        base = g.late_area.iloc[0]; target = 0.5 * base; ic50 = np.nan; prev = base
        for i in range(1, len(g)):
            cur = g.late_area.iloc[i]
            if np.isnan(cur):
                continue
            if prev >= target >= cur:
                x0, x1 = g.cu.iloc[i - 1], g.cu.iloc[i]
                ic50 = x0 + (target - prev) * (x1 - x0) / (cur - prev)
                break
            prev = cur
        ic.append({"strain_code": sc, "IC50_est": ic50})
    cr = cr.merge(pd.DataFrame(ic), on="strain_code", how="left")
    new = color.merge(cr, on="strain_code", how="outer").merge(texture, on="strain_code", how="outer")
    new.to_csv(out / "gwas_next_phenotypes_duckdb.csv", index=False)
    print(f"[pheno] new trait table: {new.shape}")

    # compare to existing old-data fam-order table
    old = pd.read_csv(REPO / "analysis/gwas/results/gwas_next_phenotypes_fam_order.csv")
    oldt = old.rename(columns={f"tex_{m}": f"tex_{m}" for m in TEX})
    j = oldt.merge(new, on="strain_code", suffixes=("_old", "_new"), how="left")
    print(f"[compare] old fam-order strains={len(old)}; present in new table={j.AUC_0_new.notna().sum() if 'AUC_0_new' in j else 'NA'}; "
          f"strains in fam missing from new table={new.strain_code.isin(old.strain_code).sum()}/{len(old)} present")
    comp = []
    for c in [c for c in old.columns if c in new.columns and c not in ("strain_code",)]:
        x, y = j[f"{c}_old"], j[f"{c}_new"]
        ok = x.notna() & y.notna()
        if ok.sum() < 5:
            comp.append({"trait": c, "n": int(ok.sum())}); continue
        comp.append({"trait": c, "n": int(ok.sum()),
                     "pearson": float(stats.pearsonr(x[ok], y[ok])[0]),
                     "spearman": float(stats.spearmanr(x[ok], y[ok])[0]),
                     "old_median": float(x[ok].median()), "new_median": float(y[ok].median())})
    comp = pd.DataFrame(comp)
    comp.to_csv(out / "trait_old_vs_new_comparison.csv", index=False)
    print(comp.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    fam_order = old.strain_code.tolist()
    al = new[new.strain_code.isin(fam_order)].set_index("strain_code").loc[[s for s in fam_order if s in set(new.strain_code)]].reset_index()
    al.to_csv(out / "gwas_next_phenotypes_fam_order_duckdb.csv", index=False)
    print(f"[pheno] aligned {len(al)}/{len(fam_order)} fam-order strains")


if __name__ == "__main__":
    main()
