#!/usr/bin/env python3
"""Idea 05 (Representation Learning): a low-dimensional atlas of colony phenotype.

Rank-PCA + varimax factor discovery over a curated 38-trait vector (shape,
intensity, CIELAB, HSV, xy, texture) at Cu=0 late window; UMAP + HDBSCAN
unsupervised atlas crossed against taxonomy/environment; pigment-factor
structure.
"""
from __future__ import annotations

import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sklearn
from sklearn.decomposition import PCA
from sklearn.preprocessing import quantile_transform, StandardScaler
import umap
import hdbscan

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import common as CM

CURATED = [
    "Shape_Area", "Shape_Perimeter", "Shape_Compactness", "Shape_Circularity",
    "Shape_Solidity", "Shape_Extent", "Shape_Eccentricity", "Shape_MeanRadius",
    "Shape_MaxRadius",
    "Intensity_MeanIntensity", "Intensity_MedianIntensity",
    "Intensity_CoefficientVarianceIntensity",
    "ColorLab_L*Median", "ColorLab_L*CoeffVar",
    "ColorLab_a*Median", "ColorLab_a*CoeffVar",
    "ColorLab_b*Median", "ColorLab_b*CoeffVar",
    "ColorLab_ChromaEstimatedMedian",
    "ColorHSV_HueMedian", "ColorHSV_HueStdDev", "ColorHSV_HueCoeffVar",
    "ColorHSV_SaturationMedian", "ColorHSV_SaturationCoeffVar",
    "ColorHSV_BrightnessMedian", "ColorHSV_BrightnessStdDev",
    "Colorxy_xMedian", "Colorxy_yMedian",
    "TextureGray_Entropy-avg-scale05", "TextureGray_Contrast-avg-scale05",
    "TextureGray_AngularSecondMoment-avg-scale05",
    "TextureGray_Correlation-avg-scale05",
    "TextureGray_InfoCorrelation1-avg-scale05",
    "TextureGray_SumEntropy-avg-scale05",
    "TextureGray_DiffEntropy-avg-scale05",
]
BLOCK = {}
for f in CURATED:
    b = "lab" if f.startswith("ColorLab") else "hsv" if f.startswith("ColorHSV") \
        else "xy" if f.startswith("Colorxy") else "texture" if f.startswith("TextureGray") \
        else "shape" if f.startswith("Shape") else "intensity"
    BLOCK[f] = b


def varimax(Phi, gamma=1.0, q=1000, tol=1e-8):
    p, k = Phi.shape
    R = np.eye(k)
    d = 0
    for _ in range(q):
        d_old = d
        Z = Phi @ R
        u, s, vh = np.linalg.svd(Phi.T @ (Z ** 3 - (gamma / p) * Z @ np.diag(np.sum(Z * Z, 0))))
        R = u @ vh
        d = np.sum(s)
        if d_old != 0 and d / d_old < 1 + tol:
            break
    return Phi @ R, R


def main() -> None:
    print("[idea05] loading extract ...")
    df = CM.read_extract()
    d = df[df.strain_id.notna() & ~df.is_control.astype(bool)
           & (df.tp_h >= 85) & (df.tp_h <= 110)].copy()
    mat = d.groupby("strain_id")[CURATED].median().astype(float)
    meta0 = CM.read_meta()[["strain_id", "species", "environment", "origin"]]
    ncol = d.groupby("strain_id").size().rename("n_colonies").reset_index()
    meta = meta0.merge(ncol, on="strain_id", how="left")
    keep = meta[meta["species"].notna()][["strain_id"]].merge(mat, left_on="strain_id",
                                                              right_index=True)
    mat = keep.drop(columns="strain_id").set_index(keep["strain_id"])
    mat = mat.dropna()
    meta = meta[meta["species"].notna()].set_index("strain_id").loc[mat.index]
    print(f"  [idea05] strains={len(mat)} traits={mat.shape[1]}")

    # rank-gaussian + standard z
    mat_rank = quantile_transform(mat, output_distribution="normal", n_quantiles=min(len(mat), 1000),
                                  random_state=42)
    Xz = StandardScaler().fit_transform(mat_rank)

    # rank-PCA
    pca = PCA(n_components=min(12, len(mat) - 1)).fit(Xz)
    evr = pd.DataFrame({"pc": np.arange(1, len(pca.explained_variance_ratio_) + 1),
                        "var_frac": pca.explained_variance_ratio_,
                        "cumvar": np.cumsum(pca.explained_variance_ratio_)})
    CM.save(evr, "idea05_pca_variance.csv")

    # varimax on top k (>=90% cumvar)
    k = int((evr.cumvar >= 0.90).idxmax()) + 1
    k = min(k, 8)
    scores = pca.components_[:k].T * np.sqrt(pca.explained_variance_[:k])
    rot, Rrot = varimax(scores)
    ldf = pd.DataFrame(rot, columns=[f"F{i+1}" for i in range(k)], index=mat.columns)
    ldf["block"] = [BLOCK[c] for c in mat.columns]
    ldf.reset_index(inplace=True)
    ldf.rename(columns={"index": "feature"}, inplace=True)
    CM.save(ldf, "idea05_varimax_loadings.csv")

    # factor scores: rotate PCA score space by the varimax rotation
    raw_scores = Xz @ pca.components_[:k].T
    fac_score = pd.DataFrame(raw_scores @ Rrot, columns=[f"F{i+1}" for i in range(k)],
                             index=mat.index)

    # block loading summary per factor
    blk = ldf.melt(id_vars=["feature", "block"], var_name="factor", value_name="loading")
    blk = blk.groupby(["factor", "block"])["loading"].apply(lambda s: s.abs().max()).unstack()
    CM.save(blk.reset_index(), "idea05_factor_block_top.csv")

    # UMAP + HDBSCAN atlas
    um = umap.UMAP(n_neighbors=15, min_dist=0.1, n_components=2, random_state=42)
    emb = um.fit_transform(Xz)
    # HDBSCAN on the UMAP projection: full-dim space is a continuum, but
    # connectivity neighborhoods in the projected manifold are separable.
    mem = hdbscan.HDBSCAN(min_cluster_size=8, min_samples=3).fit(emb)
    at = pd.DataFrame({
        "strain_id": mat.index,
        "umap_x": emb[:, 0], "umap_y": emb[:, 1],
        "hdbscan_cluster": mem.labels_, "hdbscan_prob": mem.probabilities_,
    })
    at = at.merge(meta, on="strain_id", how="left")
    CM.save(at, "idea05_atlas.csv")

    # atlas contour: cluster x species / env
    cont = (at.value_counts(["hdbscan_cluster", "species"])
            .rename("n").reset_index()
            .pivot(index="hdbscan_cluster", columns="species", values="n").fillna(0.0))
    cont["n_total"] = cont.sum(axis=1)
    CM.save(cont.reset_index(), "idea05_cluster_species_matrix.csv")

    print("  [idea05] n_clusters(hdbscan) =", len(cont) - (at.hdbscan_cluster == -1).sum(),
          f"(noise={int((at.hdbscan_cluster==-1).sum())})")

    # ---- figures ----
    plt.figure(figsize=(9, 5))
    top_cs = [c for c in cont.columns if c != "n_total"]
    per_c = cont[top_cs].div(cont["n_total"], axis=0)
    per_c.T.iloc[::-1].plot(kind="barh", stacked=True, figsize=(9, 5), cmap="tab20")
    plt.xlabel("cluster composition (fraction of species)")
    plt.title("Unsupervised phenotype clusters -> species composition")
    plt.tight_layout(); plt.savefig(CM.FIGURES / "fig05_cluster_species.png", dpi=150)
    plt.close()

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    sps = at["species"].dropna().unique()
    cmap = plt.get_cmap("tab20", len(sps))
    for i, s in enumerate(sps):
        m = at[at["species"] == s]
        axes[0].scatter(m.umap_x, m.umap_y, s=8, alpha=0.6, color=cmap(i), label=s.split()[-1])
    axes[0].set_title("UMAP colored by species"); axes[0].legend(fontsize=6, ncol=2)
    for c in sorted(x for x in at.hdbscan_cluster.unique() if x >= 0):
        m = at[at.hdbscan_cluster == c]
        axes[1].scatter(m.umap_x, m.umap_y, s=8, alpha=0.6, label=f"cluster {c}")
    axes[1].set_title("HDBSCAN clusters (noise grey)")
    axes[1].scatter(at[at.hdbscan_cluster == -1].umap_x, at[at.hdbscan_cluster == -1].umap_y,
                    s=8, alpha=0.6, color="0.7", label="noise")
    axes[1].legend(fontsize=7, ncol=3)
    plt.tight_layout(); plt.savefig(CM.FIGURES / "fig05_atlas_umap.png", dpi=150)
    plt.close()
    print("[idea05] done.")


if __name__ == "__main__":
    sys.exit(main())
