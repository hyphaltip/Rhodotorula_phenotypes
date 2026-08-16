#!/usr/bin/env python3
"""Shared data layer for the 2026-08-15 ideation session.

Extracts the per-colony-object time-course feature set from the DuckDB
v_phenotype view into a single gzipped TSV, plus strain metadata.

Rationale for the curated column set (grounded in findings up to 2026-08-15):
- Imaging "pass" is ROUND(hours_since_plate_start) -- see learning L-8
  (rounding to integer pass restores the full set of colonies imaged
  minutes apart within one pass).
- Named DuckDB parameters ':x' are NOT supported in this duckdb version;
  this script uses only literal SQL (no bind params).
- L* ~ Intensity collinearity (r>=0.99) documented; we keep both but
  remember to never regress them against each other blindly.
"""
from __future__ import annotations

import io
import pathlib
import sys

import duckdb
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
SESSION = HERE.parent  # analysis/ideas/2026-08-15-color-phenotype-space
DB = HERE.parents[3] / "db" / "rhodotorula_phenotypes.duckdb"
OUTDIR = SESSION / "data"
OUTDIR.mkdir(exist_ok=True)

SHAPE = [
    "Shape_Area", "Shape_Perimeter", "Shape_Circularity", "Shape_ConvexArea",
    "Shape_MedianRadius", "Shape_MeanRadius", "Shape_MaxRadius",
    "Shape_MinFeretDiameter", "Shape_MaxFeretDiameter", "Shape_Eccentricity",
    "Shape_Solidity", "Shape_Extent", "Shape_BboxArea",
    "Shape_MajorAxisLength", "Shape_MinorAxisLength", "Shape_Compactness",
]
INTENSITY = [
    "Intensity_IntegratedIntensity", "Intensity_MeanIntensity",
    "Intensity_MedianIntensity", "Intensity_StandardDeviationIntensity",
    "Intensity_CoefficientVarianceIntensity",
    "Intensity_InterquartileRangeIntensity", "Intensity_Density",
    "Intensity_ConvexDensity",
]
COLLAB = [c for c in (
    "ColorLab_L*Min ColorLab_L*Q1 ColorLab_L*Mean ColorLab_L*Median "
    "ColorLab_L*Q3 ColorLab_L*Max ColorLab_L*StdDev ColorLab_L*CoeffVar "
    "ColorLab_a*Min ColorLab_a*Q1 ColorLab_a*Mean ColorLab_a*Median "
    "ColorLab_a*Q3 ColorLab_a*Max ColorLab_a*StdDev ColorLab_a*CoeffVar "
    "ColorLab_b*Min ColorLab_b*Q1 ColorLab_b*Mean ColorLab_b*Median "
    "ColorLab_b*Q3 ColorLab_b*Max ColorLab_b*StdDev ColorLab_b*CoeffVar "
    "ColorLab_ChromaEstimatedMean ColorLab_ChromaEstimatedMedian").split()]
COLLHSV = []
for chan in ("Hue", "Saturation", "Brightness"):
    for stat in ("Min", "Q1", "Mean", "Median", "Q3", "Max", "StdDev", "CoeffVar"):
        COLLHSV.append(f"ColorHSV_{chan}{stat}")
COLLXY = []
for c in ("x", "y"):
    for stat in ("Min", "Q1", "Mean", "Median", "Q3", "Max", "StdDev", "CoeffVar"):
        COLLXY.append(f"Colorxy_{c}{stat}")
TEXG = [
    f"TextureGray_{name}-avg-scale05" for name in (
        "AngularSecondMoment", "Contrast", "Correlation", "HaralickVariance",
        "InverseDifferenceMoment", "SumAverage", "SumVariance", "SumEntropy",
        "Entropy", "DiffVariance", "DiffEntropy", "InfoCorrelation1",
        "InfoCorrelation2")
]

KEY = ["strain_id", "strain_code", "species", "run_number", "plate_number",
       "well_position", "object_label", "replicate_label", "is_control",
       "copper_mm", "hours_since_plate_start", "tp_h", "genus"]
ALL = KEY + SHAPE + INTENSITY + COLLAB + COLLHSV + COLLXY + TEXG


def q(col: str) -> str:
    return f'"{col}"'


def main() -> None:
    con = duckdb.connect(str(DB), read_only=True)
    cols = [q(c) for c in ALL if c not in ("copper_mm", "tp_h", "genus")]
    sel = [
        "strain_id", "strain_code", "species",
        "run_number", "plate_number", "well_position", "object_label",
        "replicate_label", "is_control",
        'CAST(factors[\'Copper concentration\'] AS DOUBLE) AS copper_mm',
        "hours_since_plate_start",
        "ROUND(hours_since_plate_start)::INT AS tp_h",
        "split_part(species,' ',1) AS genus",
    ] + [q(c) for c in SHAPE + INTENSITY + COLLAB + COLLHSV + COLLXY + TEXG]
    sql = "SELECT " + ", ".join(sel) + " FROM v_phenotype ORDER BY run_number, plate_number, strain_id, tp_h, well_position, object_label"
    df = con.execute(sql).df()
    n = len(df)
    assert n == 211800, f"expected 211800 rows from v_phenotype, got {n}"
    meta = con.execute(
        "SELECT strain_id, strain_code, strain_name, species, origin, environment FROM strain ORDER BY strain_id"
    ).df()
    con.close()

    out = OUTDIR / "db_extract.tsv.gz"
    df.to_csv(out, sep="\t", index=False, compression="gzip")
    meta_out = OUTDIR / "strain_metadata.tsv"
    meta.to_csv(meta_out, sep="\t", index=False)
    print(f"extract rows: {n}, cols: {len(ALL)}")
    print(f"rows with species: {df.species.notna().sum()} ({df.species.notna().mean():.3f})")
    print(f"control rows: {df.is_control.sum()} ({df.is_control.mean():.3f})")
    print(f"strains on extract: {df.strain_id.nunique()} (324 rows w/o strain? {df.strain_id.isna().sum()})")
    print(f"wrote {out}")


if __name__ == "__main__":
    sys.exit(main())
