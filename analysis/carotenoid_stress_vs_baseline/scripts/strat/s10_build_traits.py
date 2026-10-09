#!/usr/bin/env python3
"""Well-level colour (L*, a*, b*, chroma, hue angle) and morphology table, same wells and window as prepare_wells.py.
Largest object per well-image; absolute window [T-24, T] per metal (T = median plate span); wells need >= 2 window images.
Writes results/wells_traits.csv (main: largest object >= 2000 px, Zinc window ends 80 h) and results/wells_traits_nofilter.csv."""
import sys
import numpy as np, pandas as pd, duckdb
COLS = {"L": "ColorLab_L*GeoMedian", "a": "ColorLab_a*GeoMedian", "b": "ColorLab_b*GeoMedian", "sat": "ColorHSV_SaturationRobustMean", "val": "ColorHSV_ValueRobustMean",
        "area": "Shape_Area", "perim": "Shape_Perimeter", "circ": "Shape_Circularity", "solid": "Shape_Solidity", "ecc": "Shape_Eccentricity", "compact": "Shape_Compactness",
        "extent": "Shape_Extent", "major": "Shape_MajorAxisLength", "minor": "Shape_MinorAxisLength"}
WINDOW_H, MIN_IMG = 24.0, 2
log = lambda m: print(m, flush=True)
con = duckdb.connect("db/rhodotorula_phenotypes.duckdb", read_only=True)
sel = ", ".join(f'"{v}" as {k}' for k, v in COLS.items())
df = con.execute(f"""select Metal, run_number, plate_position, Grid_RowNum, Grid_ColNum, Concentration as conc, strain_id, capture_datetime, {sel}
                     from heavy_metal_measurement where strain_id is not null and strain_id not like 'Control%'""").df()
log(f"rows (named strains): {len(df):,}"); assert len(df) == 502051
df = df[df.Metal.isin(["Chromium", "Copper", "Lead"])].copy(); log(f"metals kept (Iron and Zinc held out, D-56): {len(df):,}")
pk = ["Metal", "run_number", "plate_position"]; wk = pk + ["Grid_RowNum", "Grid_ColNum"]
df["h"] = (df.capture_datetime - df.groupby(pk).capture_datetime.transform("min")).dt.total_seconds() / 3600; assert 90 < df.h.max() < 130
df = df.sort_values("area", ascending=False).drop_duplicates(wk + ["capture_datetime"])
T = df.groupby(pk).h.max().groupby("Metal").median(); T["Zinc"] = 80.0; df["T"] = df.Metal.map(T)   # Zinc plates stop at 83.7-89.7 h
def build(d, name):
    win = d[(d.h >= d["T"] - WINDOW_H) & (d.h <= d["T"])]
    agg = {k: (k, "median") for k in COLS}
    w = win.groupby(wk).agg(strain_id=("strain_id", "first"), conc=("conc", "first"), n_img=("capture_datetime", "nunique"), **agg).reset_index()
    w = w[w.n_img >= MIN_IMG].copy()
    w["lnA"] = np.log(w.area); w["chroma"] = np.hypot(w.a, w.b); w["hue_deg"] = np.degrees(np.arctan2(w.b, w.a)); w["aspect"] = w.major / w.minor
    w["plate_id"] = w.Metal + "_" + w.run_number.astype(str) + "_" + w.plate_position.astype(str); w["dose_s"] = w.conc / w.groupby("Metal").conc.transform("max")
    w["lnA_c"] = w.lnA - w.Metal.map(w[w.conc == 0].groupby("Metal").lnA.mean())
    bad = w[["a", "b", "circ", "solid", "ecc", "compact", "extent", "aspect"]].isna().sum(); log(f"{name}: wells {len(w):,}; NaN counts {bad[bad > 0].to_dict()}")
    w.to_csv(f"analysis/carotenoid_stress_vs_baseline/results/{name}.csv", index=False)
build(df[df.area >= 2000], "wells_traits"); build(df, "wells_traits_nofilter")
