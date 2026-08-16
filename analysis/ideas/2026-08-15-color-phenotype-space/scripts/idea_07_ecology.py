#!/usr/bin/env python3
"""Idea 07 (Trait-Ecology Spectrum): do phenotype traits segregate by isolation
habitat, and does the association survive species-stratified permutation?

Strain-level traits (Cu=0 late window): pigmentation (a*Median, chroma, L*Median),
colony size (area), heterogeneity (texture entropy, b*CoeffVar).
Environments: marsh/tidalflat, soil, plant, food, other.
Two nulls: raw label-shuffle permutation, and within-species block shuffle.
"""
from __future__ import annotations

import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import common as CM

TRAITS = {"a_median": "ColorLab_a*Median", "chroma": "ColorLab_ChromaEstimatedMedian",
          "L_median": "ColorLab_L*Median", "area": "Shape_Area",
          "tex_entropy": "TextureGray_Entropy-avg-scale05",
          "b_cv": "ColorLab_b*CoeffVar"}


def oneway_f(y: np.ndarray, g: np.ndarray) -> float:
    groups = [y[g == u] for u in np.unique(g) if np.sum(g == u) > 1]
    if len(groups) < 2:
        return 0.0
    f, _ = stats.f_oneway(*groups)
    return float(f)


def eta_sq(y: np.ndarray, g: np.ndarray) -> float:
    gm = y.mean()
    ssb = sum(len(y[g == u]) * (y[g == u].mean() - gm) ** 2 for u in np.unique(g))
    return ssb / ((y - gm) ** 2).sum()


def main() -> None:
    print("[idea07] loading extract ...")
    df = CM.read_extract()
    d = df[df.strain_id.notna() & ~df.is_control.astype(bool)
           & (df.tp_h >= 85) & (df.tp_h <= 110)].copy()
    meta = CM.read_meta()[["strain_id", "species", "environment"]].rename(
        columns={"species": "species_meta", "environment": "environment_meta"})
    d = d.merge(meta, on="strain_id", how="left")
    d["species"] = d["species_meta"].fillna(d.get("species"))
    env = d["environment_meta"].fillna("other").copy()
    d["env_group"] = np.where(env.isin(["marsh_tidalflat", "soil", "plant", "food"]),
                              env, "other")
    # collapse tiny cats to other but keep big 4
    rows = []
    rng = np.random.default_rng(7)
    for tname, feats in TRAITS.items():
        sub = d[["strain_id", "env_group", "species", feats]].dropna()
        strain = (sub.groupby("strain_id").agg(
            val=(feats, "median"), env_group=("env_group", "first"),
            species=("species", "first")).reset_index())
        g = pd.Categorical(strain.env_group).codes
        y = strain.val.values.astype(float)
        groups = np.unique(g)
        if len(groups) < 2:
            continue
        fobs = oneway_f(y, g)
        p_raw = (1 + sum(oneway_f(y, rng.permutation(g)) >= fobs for _ in range(2000))) / 2001
        # species-stratified permutation
        nz = 0
        for _ in range(2000):
            perm = g.copy()
            for sp in strain.species.dropna().unique():
                ix = np.where(strain.species == sp)[0]
                perm[ix] = g[rng.permutation(ix)]
            if oneway_f(y, perm) >= fobs:
                nz += 1
        p_sp = (1 + nz) / 2001
        rows.append({"trait": tname, "feature": feats, "n_strains": len(strain),
                     "F_raw": fobs, "p_raw": p_raw, "p_species_stratified": p_sp,
                     "eta2": eta_sq(y, g)})
    out = pd.DataFrame(rows).sort_values("p_species_stratified")
    CM.save(out, "idea07_trait_environment.csv")
    print("  [idea07] trait-environment association (F, p_raw, p_species_stratified):")
    print(out[["trait", "F_raw", "p_raw", "p_species_stratified", "eta2"]].to_string(index=False))

    # z-scored trait means by environment (typed strains only)
    envmeans = []
    for envg in ["marsh_tidalflat", "soil", "plant", "food", "other"]:
        m = d[d.env_group == envg]
        mm = m.groupby("strain_id")[list(TRAITS.values())].median()
        envmeans.append({"environment": envg, **mm.mean().to_dict()})
    em = pd.DataFrame(envmeans)
    zv = em[list(TRAITS.values())].apply(lambda s: (s - s.mean()) / s.std(), axis=0)
    CM.save(em, "idea07_env_trait_means.csv")
    CM.save(zv.assign(environment=em.environment), "idea07_env_trait_z.csv")

    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    zvT = zv.copy()
    zvT.insert(0, "environment", em.environment)
    zvTm = zvT.melt(id_vars="environment", var_name="trait", value_name="z")
    zz = zvTm.pivot(index="trait", columns="environment", values="z")
    im = ax[0].imshow(zz.values, cmap="RdBu_r", vmin=-1.6, vmax=1.6, aspect="auto")
    ax[0].set_xticks(range(len(zz.columns))); ax[0].set_xticklabels(zz.columns, rotation=30, fontsize=8)
    ax[0].set_yticks(range(len(zz.index))); ax[0].set_yticklabels(zz.index, fontsize=8)
    ax[0].set_title("trait z-score by environment (strain means, Cu=0)")
    fig.colorbar(im, ax=ax[0])
    mh = d[d.env_group == "marsh_tidalflat"]
    ot = d[d.env_group != "marsh_tidalflat"]
    for s, col, lab in [(mh, "#2e7d32", "marsh/tidalflat"), (ot, "#8e8e8e", "other env")]:
        mm = s.groupby("strain_id")["ColorLab_a*Median"].median().dropna()
        ax[1].hist(mm, bins=20, alpha=0.55, color=col, label=f"{lab} (n={len(mm)})")
    ax[1].set_xlabel("a*Median (redness)")
    ax[1].set_title("marsh/tidalflat vs other environments"); ax[1].legend()
    plt.tight_layout(); plt.savefig(CM.FIGURES / "fig07_trait_environment.png", dpi=150)
    plt.close()
    print("[idea07] done.")


if __name__ == "__main__":
    sys.exit(main())
