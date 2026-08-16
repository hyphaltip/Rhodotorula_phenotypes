#!/usr/bin/env python3
"""Idea 03 (Information Theorist): what actually carries signal in ~200-col space.

(3a) Redundancy structure + greedy forward conditional-MI selection toward a
     minimal, species-discriminating feature set; effective rank of the space.
(3b) Does Hue/Texture carry species information that L* (and a*,b*) does not?
     Conditional mutual information I(F; Species | L) etc.
"""
from __future__ import annotations

import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import common as C

FEATS = {
    "shape": ["Shape_Area", "Shape_Compactness", "Shape_Circularity", "Shape_Solidity",
              "Shape_Extent", "Shape_Eccentricity", "Shape_MajorAxisLength",
              "Shape_MeanRadius", "Shape_MaxRadius", "Shape_Perimeter"],
    "intensity": ["Intensity_MeanIntensity", "Intensity_MedianIntensity",
                  "Intensity_CoefficientVarianceIntensity",
                  "Intensity_StandardDeviationIntensity", "Intensity_Density",
                  "Intensity_InterquartileRangeIntensity"],
    "lab": ["ColorLab_L*Median", "ColorLab_a*Median", "ColorLab_b*Median",
            "ColorLab_CHROMA", "ColorLab_a*CoeffVar", "ColorLab_b*CoeffVar",
            "ColorLab_L*CoeffVar"],
    "hsv": ["ColorHSV_SaturationMedian", "ColorHSV_BrightnessMedian",
            "ColorHSV_HueMedian", "ColorHSV_SaturationCoeffVar"],
    "xy": ["Colorxy_xMedian", "Colorxy_yMedian"],
    "texture": ["TextureGray_Contrast-avg-scale05", "TextureGray_Entropy-avg-scale05",
                "TextureGray_AngularSecondMoment-avg-scale05",
                "TextureGray_Correlation-avg-scale05", "TextureGray_DiffVariance-avg-scale05",
                "TextureGray_HaralickVariance-avg-scale05",
                "TextureGray_SumEntropy-avg-scale05", "TextureGray_InfoCorrelation1-avg-scale05"],
}


def build_strain_matrix(df: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    d = df[df.strain_id.notna() & ~df.is_control.astype(bool)
           & (df.tp_h >= 85) & (df.tp_h <= 110)].copy()
    d["ColorLab_CHROMA"] = d["ColorLab_ChromaEstimatedMedian"]
    feats = [f for block in FEATS.values() for f in block]
    rows = []
    for sid, g in d.groupby("strain_id"):
        if len(g) < 5:
            continue
        r = {"strain_id": sid, "species": g.species.dropna().iloc[0] if g.species.notna().any() else np.nan,
             "n_colonies": len(g)}
        for f in feats:
            r[f] = float(np.nanmedian(g[f]))
        rows.append(r)
    m = pd.DataFrame(rows)
    m = m.merge(meta[["strain_id", "origin", "environment"]], on="strain_id", how="left")
    return m


def mi_counts(x_disc: np.ndarray, y_disc: np.ndarray, alpha: float = 0.5) -> float:
    """Empirical mutual information (bits) between two discrete vectors."""
    n = len(x_disc)
    if n == 0:
        return float("nan")
    kx, ky = np.max(x_disc) + 1, np.max(y_disc) + 1
    px = np.bincount(x_disc.astype(int), minlength=kx) + alpha
    py = np.bincount(y_disc.astype(int), minlength=ky) + alpha
    pxy = np.zeros((kx, ky))
    for (xi, yi) in zip(x_disc, y_disc):
        pxy[int(xi), int(yi)] += 1
    pxy += alpha
    pxy /= pxy.sum(); px = px / px.sum(); py = py / py.sum()
    tot = 0.0
    for i in range(kx):
        for j in range(ky):
            if pxy[i, j] > 0:
                tot += pxy[i, j] * np.log2(pxy[i, j] / (px[i] * py[j]))
    return float(tot)


def discretize(x: pd.Series, q: int) -> np.ndarray:
    x = x.astype(float)
    edges = np.nanquantile(x, np.linspace(0, 1, q + 1))
    edges[0], edges[-1] = -np.inf, np.inf
    edges = np.unique(edges)
    return np.digitize(x, edges[1:-1])


def cond_mi(x: pd.Series, y_disc: np.ndarray, cond: pd.Series, q_all: int = 4) -> float:
    """I(x; y | cond) via sum_c p(c) I(x_c; y_c) with x discretized within cond bins."""
    xd = discretize(x, q_all)
    cd = discretize(cond, q_all)
    tot = 0.0
    for c in np.unique(cd):
        m = cd == c
        if m.sum() < 10:
            continue
        w = m.sum() / len(cd)
        tot += w * mi_counts(xd[m], y_disc[m])
    return float(tot)


def forward_select(m: pd.DataFrame, species: np.ndarray, feats: list, max_k: int = 4, q: int = 4) -> list:
    from itertools import product
    disc = {f: discretize(m[f], q) for f in feats}
    sel, remaining = [], feats[:]
    for _ in range(max_k):
        best, best_i = None, -1
        for f in remaining:
            if len(sel) == 0:
                mi = mi_counts(disc[f], species)
            elif len(sel) == 1:
                joint = disc[sel[0]] * q + disc[f]
                mi = mi_counts(joint, species)
            else:
                # greedy approximation using pairwise-average conditional
                mi = np.mean([cond_mi(m[f], species, m[s], 3) for s in sel])
            if best is None or mi > best_i:
                best, best_i = f, mi
        sel.append(best)
        remaining.remove(best)
        print(f"    step {len(sel)}: +{best}  MI-combo={best_i:.3f} bits")
    return sel


def main() -> None:
    print("[idea03] loading extract ...")
    df = C.read_extract()
    meta = C.read_meta()
    print("[idea03] build strain-level matrix (Cu=0, late window) ...")
    m = build_strain_matrix(df, meta)
    m = m[m.species.notna()].copy()
    # species groups with >= 6 strains
    vc = m.species.value_counts()
    keep = vc[vc >= 6].index
    m = m[m.species.isin(keep)].copy()
    print(f"  strains={len(m)} species-groups={m.species.nunique()} "
          f"(typed strains with >=6 isolates)")
    feats = [f for block in FEATS.values() for f in block]
    s_enc, s_lab = pd.factorize(m.species)
    n_species = len(s_lab)
    H = np.log2(n_species)
    print(f"  H(species) = {H:.3f} bits")

    # per-feature MI and conditional-MI
    print("[idea03] per-feature MI (plain, |L, |a*,b*) ...")
    rows = []
    for f in feats:
        mi = mi_counts(discretize(m[f], 6), s_enc)
        m_l = cond_mi(m[f], s_enc, m["ColorLab_L*Median"], 4)
        m_ab = cond_mi(m[f], s_enc, m["ColorLab_a*Median"] + 256 * m["ColorLab_b*Median"], 4)
        rows.append({"feature": f, "block": _block(f), "MI_species_bits": mi,
                     "MI_species_fracH": mi / H if H > 0 else np.nan,
                     "CMI|L_bits": m_l, "CMI|L_a_b_bits": m_ab,
                     "info_beyond_L_ab": mi - m_ab})
    info = pd.DataFrame(rows).sort_values("MI_species_bits", ascending=False)
    C.save(info, "idea03_feature_species_info.csv")

    # redundancy / effective rank
    print("[idea03] redundancy structure ...")
    X = m[feats].apply(lambda s: s.rank(pct=True)).to_numpy()  # Spearman approx via ranks
    A = np.corrcoef(X.T)
    eig = np.linalg.eigvalsh(A)
    eig = eig[::-1]
    cum = np.cumsum(eig) / np.sum(eig)
    eff95 = int(np.searchsorted(cum, 0.95) + 1)
    rr = []
    for i, a in enumerate(feats):
        for j in range(i + 1, len(feats)):
            rr.append((abs(A[i, j]), feats[i], feats[j]))
    top = sorted(rr, reverse=True)[:10]
    red = pd.DataFrame({"effective_rank_95var": eff95, "n_features": len(feats),
                        "mean_abs_corr": float(np.mean(np.abs(A[np.triu_indices(len(feats), 1)])))},
                       index=[0])
    C.save(red, "idea03_redundancy_summary.csv")
    C.save(pd.DataFrame(top, columns=["abs_corr", "f1", "f2"]), "idea03_top_correlated_pairs.csv")

    # greedy forward selection toward minimal species-discriminating set
    print("[idea03] greedy forward conditional-MI selection ...")
    sel = forward_select(m, s_enc, feats, max_k=4)
    C.save(pd.DataFrame({"selected_step": range(1, len(sel) + 1), "feature": sel}),
           "idea03_forward_selection.csv")

    # ---------- figures ----------
    plt.figure(figsize=(10, 6))
    topf = info.head(25)
    colors = {"lab": "#d62728", "hsv": "#ff7f0e", "xy": "#9467bd", "texture": "#2ca02c",
              "shape": "#1f77b4", "intensity": "#e377c2"}
    bars = plt.bar(np.arange(len(topf)), topf.MI_species_bits,
                   color=[colors[b] for b in topf.block])
    plt.xticks(np.arange(len(topf)), [f.split()[1] if " " in f else f.split("_")[0] for f in topf.feature],
               rotation=60, fontsize=6)
    plt.ylabel("MI(feature; species) [bits]")
    plt.title("Species information carried by each phenotype feature (Cu=0, late)")
    plt.tight_layout(); plt.savefig(C.FIGURES / "fig03_mi_top.png", dpi=150)
    plt.close()

    # heatmap correlation matrix (clustered)
    Z = linkage(1 - A, method="average")
    order = list(leaves_list(Z))
    Ao = A[np.ix_(order, order)]
    lo = [feats[i] for i in order]
    fig, ax = plt.subplots(figsize=(8.5, 7.5))
    im = ax.imshow(Ao, vmin=-1, vmax=1, cmap="RdBu_r", aspect="auto")
    ax.set_xticks(range(len(lo))); ax.set_xticklabels(lo, rotation=90, fontsize=4.5)
    ax.set_yticks(range(len(lo))); ax.set_yticklabels(lo, fontsize=4.5)
    fig.colorbar(im, label="Spearman r")
    ax.set_title("Strain-level feature correlation matrix (clustered)")
    plt.tight_layout(); plt.savefig(C.FIGURES / "fig03_redundancy_heatmap.png", dpi=150)
    plt.close()
    print("[idea03] done.")


def _block(f):
    for b, feats in FEATS.items():
        if f in feats:
            return b
    return "other"


if __name__ == "__main__":
    sys.exit(main())
