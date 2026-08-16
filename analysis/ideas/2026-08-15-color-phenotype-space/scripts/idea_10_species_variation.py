#!/usr/bin/env python3
"""Idea 10 (Within-Species Variation + Scale of Variation).

Build a per-strain trait table (the 5 traits idea 09 tested on the tree),
then summarise the SCALE of variation:
  1. Box-and-whisker plots per species (species with n>=3 strains) for each trait
     -> visual within-species spread vs between-species location shift.
  2. Variance decomposition per trait: total variance of strain values split into
     between-species variance (variance of species means) and within-species
     variance (mean of within-species variances) -> what fraction of the total
     variation is inside species, i.e. among strains of the same species?
  3. Per-species summary stats (mean, SD, CV%, range, IQR) to quantify the
     spread within each species.
"""
from __future__ import annotations

import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import common as CM

TRAITS = [
    ("slope_logchroma_per_mM", "Cu sensitivity slope", "slope per mM"),
    ("intercept_logchroma", "Baseline chroma (Cu=0)", "log10(chroma)"),
    ("l10med_fixed", "Colony size", "log10(area px)"),
    ("partial_slope_sd_cu", "Within-strain heterogeneity widening", "sd slope per mM"),
    ("pace_loglog", "Pigment pace", "dlogC/dlogA"),
]
MIN_N_SPECIES = 3


def build_trait_table() -> pd.DataFrame:
    rn = pd.read_csv(CM.RESULTS / "idea04_reaction_norms.csv")
    ws = pd.read_csv(CM.RESULTS / "idea01b_within_strain.csv")
    pp = pd.read_csv(CM.RESULTS / "idea08_pigment_pace.csv")

    pace = (pp.groupby(["strain_id", "species"], as_index=False)["pace_loglog"].median())

    t = (rn[["strain_id", "species", "slope_logchroma_per_mM", "intercept_logchroma"]]
           .merge(ws[["strain_id", "l10med_fixed", "partial_slope_sd_cu"]], on="strain_id", how="outer")
           .merge(pace, on=["strain_id", "species"], how="outer"))
    return t


def variance_decomposition(t: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for col, _lab, _un in TRAITS:
        d = t[["species", col]].dropna()
        if len(d) < 10:
            continue
        nsp = d.groupby("species")[col].agg(["count", "mean"])
        nsp = nsp[nsp["count"] >= MIN_N_SPECIES]
        grand = d[col].mean()
        n = len(d)
        # exact SS decomposition: SS_total = SS_between + SS_within
        ss_within = ((d[col] - d.groupby("species")[col].transform("mean")) ** 2).sum()
        ss_between = (nsp["count"] * (nsp["mean"] - grand) ** 2).sum()
        ss_total = ss_within + ss_between
        dfg, dfe = len(nsp) - 1, n - len(nsp)
        ms_between = ss_between / dfg if dfg > 0 else np.nan
        ms_within = ss_within / dfe if dfe > 0 else np.nan
        rows.append({
            "trait": col,
            "n_strains": n,
            "n_species": len(nsp),
            "ss_total": ss_total,
            "ss_between": ss_between,
            "ss_within": ss_within,
            "frac_between": ss_between / ss_total,
            "frac_within": ss_within / ss_total,
            "eta2_between": ss_between / ss_total,
            "F": ms_between / ms_within if ms_within > 0 else np.nan,
            "sd_total": np.sqrt(ss_total / (n - 1)),
            "sd_within_mean": np.sqrt(ms_within),
            "within_sd_as_pct_of_total": np.sqrt(ms_within) / np.sqrt(ss_total / (n - 1)) * 100,
        })
    return pd.DataFrame(rows)


def species_summary(t: pd.DataFrame, decomp: pd.DataFrame) -> pd.DataFrame:
    keep = decomp.trait.unique()
    rows = []
    for col, _lab, _un in TRAITS:
        if col not in decomp.trait.values:
            continue
        for sp, g in t.groupby("species"):
            v = g[col].dropna()
            if len(v) < MIN_N_SPECIES:
                continue
            q1, q3 = v.quantile([0.25, 0.75])
            rows.append({
                "trait": col, "species": sp, "n": len(v),
                "mean": v.mean(), "median": v.median(), "sd": v.std(ddof=1),
                "cv_pct": v.std(ddof=1) / abs(v.mean()) * 100 if v.mean() else np.nan,
                "iqr": q3 - q1, "range": v.max() - v.min(),
            })
    return pd.DataFrame(rows)


def make_boxplots(t: pd.DataFrame, decomp: pd.DataFrame) -> None:
    traitlevels = {c: (lab, un) for c, lab, un in TRAITS}
    cmap = plt.get_cmap("tab20")
    for col, (lab, un) in traitlevels.items():
        d = t[["species", col]].dropna()
        nsp = d.groupby("species")[col].count()
        d = d[d["species"].isin(nsp[nsp >= MIN_N_SPECIES].index)].copy()
        order = d.groupby("species")[col].median().sort_values().index
        d["sp"] = pd.Categorical(d["species"], categories=order, ordered=True)
        # fix colors by rank of median so hue is consistent across panels
        rank = {s: i for i, s in enumerate(order)}
        d["spi"] = d["species"].map(rank)
        d["col"] = d["spi"].map(lambda i: cmap(i / max(len(order) - 1, 1)))

        fig, ax = plt.subplots(figsize=(max(6.5, 0.55 * len(order) + 1.5), 4.6))
        for (s, g), c in zip(d.groupby("species", sort=False), d.sort_values("spi").col.unique()):
            bp = ax.boxplot(g[col].values, positions=[g.spi.iloc[0]],
                            widths=0.62, patch_artist=True,
                            medianprops=dict(color="black", lw=1.2),
                            flierprops=dict(markersize=3))
            bp["boxes"][0].set_facecolor(c)
            bp["boxes"][0].set_alpha(0.55)
            n = len(g)
            ax.text(g.spi.iloc[0], g[col].min() - 0.03 * np.ptp(d[col]), f"n={n}",
                    ha="center", va="top", fontsize=7)
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels([s.replace("Rhodotorula ", "").replace("sp. clade ", "clade ")
                            for s in order], rotation=55, ha="right", fontsize=8)
        row = decomp.loc[decomp.trait.eq(col), ["frac_between", "within_sd_as_pct_of_total"]].iloc[0]
        ax.set_xlabel("")
        ax.set_ylabel(f"{lab}\n[{un}]")
        ax.set_title(f"{lab} by species  |  between-species {row.frac_between*100:.0f}% of "
                     f"variance, within-species SD ≈ {row.within_sd_as_pct_of_total:.0f}% of total SD",
                     fontsize=9)
        ax.grid(axis="y", ls=":", alpha=0.4)
        fig.tight_layout()
        out = CM.FIGURES / f"fig10_boxplot_{col}.png"
        fig.savefig(out, dpi=180)
        plt.close(fig)
        print(f"  wrote {out.name}")


def main() -> None:
    print("[idea10] building per-strain trait table ...")
    t = build_trait_table()
    print(f"  {len(t)} strains x {len(TRAITS)} traits")

    decomp = variance_decomposition(t)
    CM.save(decomp, "idea10_variance_decomposition.csv")

    summ = species_summary(t, decomp)
    CM.save(summ, "idea10_species_variation.csv")

    make_boxplots(t, decomp)

    print("\n=== Variance decomposition (scale of variation) ===")
    print(decomp.round(4).to_string(index=False))
    print("\n=== Per-species variation (top by within-species SD) ===")
    top = (summ.sort_values("sd", ascending=False)
               .head(15)[["trait", "species", "n", "mean", "sd", "cv_pct", "iqr", "range"]])
    print(top.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
