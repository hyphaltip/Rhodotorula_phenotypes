#!/usr/bin/env python3
"""Shared helpers for the 2026-08-15 ideation analyses."""
from __future__ import annotations

import pathlib

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
SESSION = HERE.parent
DATA = SESSION / "data"
RESULTS = SESSION / "results"
FIGURES = SESSION / "figures"
RESULTS.mkdir(exist_ok=True)
FIGURES.mkdir(exist_ok=True)

PARQUET = DATA / "db_extract.parquet"
META = DATA / "strain_metadata.tsv"

COLLAB = [
    "ColorLab_L*Median", "ColorLab_a*Median", "ColorLab_b*Median",
    "ColorLab_L*Q1", "ColorLab_L*Q3", "ColorLab_L*StdDev",
    "ColorLab_L*CoeffVar", "ColorLab_a*CoeffVar", "ColorLab_b*CoeffVar",
    "ColorLab_ChromaEstimatedMean", "ColorLab_ChromaEstimatedMedian",
]
COLLHSV = [f"ColorHSV_{c}{s}" for c in ("Saturation", "Brightness", "Hue")
           for s in ("Median", "StdDev", "CoeffVar")]
_tex = ["Contrast", "Entropy", "AngularSecondMoment", "Correlation",
        "DiffVariance", "DiffEntropy", "HaralickVariance", "InverseDifferenceMoment"]
TEXG = [f"TextureGray_{n}-avg-scale05" for n in _tex]
SHAPE = ["Shape_Area", "Shape_Perimeter", "Shape_Circularity", "Shape_ConvexArea",
         "Shape_MedianRadius", "Shape_MaxRadius", "Shape_Eccentricity",
         "Shape_Solidity", "Shape_Extent", "Shape_MajorAxisLength",
         "Shape_MinorAxisLength", "Shape_Compactness"]


def read_extract() -> pd.DataFrame:
    df = pd.read_parquet(PARQUET)
    return df


def read_meta() -> pd.DataFrame:
    return pd.read_csv(META, sep="\t")


def save(df: pd.DataFrame, name: str) -> None:
    out = RESULTS / name
    df.to_csv(out, index=False)
    print(f"  wrote {out.name} [{df.shape[0]} rows]")


def boot_ci(x: np.ndarray, stat_fn, n_boot: int = 999, seed: int = 123) -> tuple[float, float, float]:
    """Bootstrap 95% CI (percentile) for stat_fn over sample x."""
    rng = np.random.default_rng(seed)
    est = stat_fn(x)
    if len(x) == 0:
        return np.nan, np.nan, np.nan
    bs = np.array([stat_fn(x[rng.integers(0, len(x), len(x))]) for _ in range(n_boot)])
    lo, hi = np.nanpercentile(bs, 2.5), np.nanpercentile(bs, 97.5)
    return float(est), float(lo), float(hi)


def ci_str(est_lo_hi: tuple) -> str:
    e, lo, hi = est_lo_hi
    if np.isnan(e):
        return "n/a"
    return f"{e:.3f} [{lo:.3f}, {hi:.3f}]"


def hue_deg_from_ab(a: pd.Series, b: pd.Series) -> pd.Series:
    """CIELAB hue angle (0-360) from per-colony a*, b* medians; NaN where chroma ~ 0."""
    h = (np.degrees(np.arctan2(b, a))) % 360.0
    chroma = np.hypot(a, b)
    return pd.Series(np.where(chroma < 1e-6, np.nan, h), index=a.index)


def circular_stats(x_deg: np.ndarray, w: np.ndarray | None = None) -> dict:
    """Chroma-weighted (or uniform) resultant vector stats on angles (deg)."""
    t = np.radians(np.asarray(x_deg, dtype=float))
    m = np.isfinite(t) & (np.asarray(x_deg) == np.asarray(x_deg))
    t = t[m]
    w = np.asarray(w, dtype=float)[m] if w is not None else np.ones_like(t)
    w = w / np.sum(w)
    sx, sy = np.sum(w * np.cos(t)), np.sum(w * np.sin(t))
    R = np.hypot(sx, sy)
    mu = np.degrees(np.arctan2(sy, sx)) % 360.0
    return {"mean_deg": float(mu), "concentration_R": float(R), "n": int(m.sum())}
