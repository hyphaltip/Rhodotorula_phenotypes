#!/usr/bin/env python3
"""Build one row per replicate well from heavy_metal_measurement.

Well = (Metal, run_number, plate_position, Grid_RowNum, Grid_ColNum). Each plate carries one concentration.
Steps (row counts are logged):
  1. drop NULL strain_id and Control-N rows
  2. per well x image keep the LARGEST object (other objects in the same well are treated as fragments)
  3. late window = ABSOLUTE hours since plate start: [T-WINDOW_H, T], T = median plate imaging span of that metal
     (plates differ in span; a plate-relative window would compare different growth stages)
  4. per well: median a*, median ln(area) over window images (need >= MIN_IMG distinct images)
"""
import argparse, sys
from pathlib import Path
import numpy as np, pandas as pd, duckdb

A = "ColorLab_a*GeoMedian"
WINDOW_H, MIN_IMG = 24.0, 2
log = lambda m: print(m, flush=True)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", default="db/rhodotorula_phenotypes.duckdb")
    ap.add_argument("--out", default="analysis/carotenoid_stress_vs_baseline/results")
    ap.add_argument("--min-area", type=float, default=0.0, help="drop well-images whose chosen object has Shape_Area < this (pixels); handles tiny colonies whose a* is background-dominated")
    ap.add_argument("--all-objects", action="store_true", help="sensitivity: median over all objects instead of the largest")
    a = ap.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect(a.db, read_only=True)
    df = con.execute(f'''select Metal, run_number, plate_position, Grid_RowNum, Grid_ColNum, Concentration, strain_id,
                         capture_datetime, Shape_Area, "{A}" as astar from heavy_metal_measurement''').df()
    log(f"1. loaded rows: {len(df):,}"); assert len(df) == 517371
    assert df[["Grid_RowNum", "Grid_ColNum", "plate_position", "run_number"]].notna().all().all(), "NULL well keys"
    df = df[df.strain_id.notna() & ~df.strain_id.str.startswith("Control")].copy()
    log(f"   after dropping NULL strain_id and Control-N: {len(df):,}")
    assert df.astar.notna().all() and (df.Shape_Area > 0).all()
    pk = ["Metal", "run_number", "plate_position"]; wk = pk + ["Grid_RowNum", "Grid_ColNum"]
    # a well must belong to exactly one strain and one dose
    chk = df.groupby(wk).agg(ns=("strain_id", "nunique"), nd=("Concentration", "nunique"))
    assert (chk.ns == 1).all() and (chk.nd == 1).all(), f"wells with >1 strain or dose: {(chk.ns>1).sum()}, {(chk.nd>1).sum()}"
    log(f"   wells: {len(chk):,}")
    df["lnA"] = np.log(df.Shape_Area)
    df["h"] = (df.capture_datetime - df.groupby(pk).capture_datetime.transform("min")).dt.total_seconds() / 3600
    df["hmax"] = df.groupby(pk).h.transform("max")
    assert 90 < df.h.max() < 130, f"time axis wrong: {df.h.max()}"
    iw = wk + ["capture_datetime"]
    df["n_obj"] = df.groupby(iw).astar.transform("size")
    if not a.all_objects:
        df = df.sort_values("Shape_Area", ascending=False).drop_duplicates(iw).sort_index()
    log(f"2. {'all objects' if a.all_objects else 'largest object per well-image'}: {len(df):,} rows; mean objects/well-image before: see n_obj")
    T = df.groupby(pk).h.max().groupby("Metal").median()
    log("   window end T per metal (median plate span, h): " + ", ".join(f"{m}={t:.1f}" for m, t in T.items()))
    df["T"] = df.Metal.map(T)
    if a.min_area > 0:
        n_before = len(df); df = df[df.Shape_Area >= a.min_area]
        log(f"   min-area filter {a.min_area:g} px: kept {len(df):,} of {n_before:,} well-images ({len(df)/n_before:.1%})")
    win = df[(df.h >= df["T"] - WINDOW_H) & (df.h <= df["T"])]
    log(f"3. late window [T-{WINDOW_H:.0f}, T] h since plate start: {len(win):,} of {len(df):,} rows")
    assert 0.05 * len(df) < len(win) < 0.5 * len(df)
    n_all = df.groupby(wk).size().reset_index(name="x").merge(df.groupby(wk).agg(conc=("Concentration", "first")).reset_index())
    w = (win.groupby(wk).agg(strain_id=("strain_id", "first"), conc=("Concentration", "first"), n_img=("capture_datetime", "nunique"),
                              a=("astar", "median"), lnA=("lnA", "median"), mean_n_obj=("n_obj", "mean")).reset_index())
    log(f"4. wells with window data: {len(w):,} of {len(n_all):,} wells with any data; "
        f"dropped for < {MIN_IMG} window images: {(w.n_img < MIN_IMG).sum():,}")
    kept = set(map(tuple, w[w.n_img >= MIN_IMG][wk].values))
    n_all["kept"] = [tuple(r) in kept for r in n_all[wk].values]
    lost = n_all.groupby(["Metal", "conc"]).agg(wells=("kept", "size"), kept=("kept", "sum")).reset_index()
    lost["lost"] = lost.wells - lost.kept
    lost.to_csv(out / "wells_lost_by_dose.csv", index=False)
    log("   wells lost by metal x dose (survivor-selection check; wells that never produced objects are not in the data at all):")
    log(lost[lost.lost > 0].to_string(index=False) if (lost.lost > 0).any() else "   none lost")
    w = w[w.n_img >= MIN_IMG].copy()
    w["plate_id"] = w.Metal + "_" + w.run_number.astype(str) + "_" + w.plate_position.astype(str)
    w["dose_s"] = w.conc / w.groupby("Metal").conc.transform("max")
    ctrl_mean = w[w.conc == 0].groupby("Metal").lnA.mean()
    w["lnA_c"] = w.lnA - w.Metal.map(ctrl_mean)
    for m, g in w.groupby("Metal"):
        r = g.groupby(["strain_id", "conc"]).size()
        log(f"   {m}: {len(g):,} wells, {g.strain_id.nunique()} strains, replicate wells per strain x dose: "
            f"median {r.median():.0f}, min {r.min()}, max {r.max()}; dose levels {sorted(g.conc.unique())}")
    fn = out / ("wells_allobj.csv" if a.all_objects else ("wells.csv" if a.min_area == 0 else f"wells_min{a.min_area:g}.csv"))
    w.to_csv(fn, index=False); log(f"wrote {fn}")

if __name__ == "__main__":
    sys.exit(main())
