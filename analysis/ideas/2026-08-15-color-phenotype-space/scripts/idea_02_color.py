#!/usr/bin/env python3
"""Idea 02 (Color/Imaging Scientist): colorimetry done right.

(2a) Chroma-weighted hue space: species pigment identity (CIELAB hue angle) vs
     pigment amount (chroma) vs tissue density (L*); stress-induced intra-colony
     heterogeneity (a*CoeffVar); pigment morph clustering on the circle.
(2b) Pigment appearance time: per-run color calibration + onset/darkening times
     across the 0-117 h time-course.
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
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from scipy.stats import spearmanr

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import common as C

LATE = (80, 110)
RNG = np.random.default_rng(7)


def prep(df: pd.DataFrame) -> pd.DataFrame:
    df = df[df.strain_id.notna() & ~df.is_control.astype(bool)].copy()
    df["chroma"] = df["ColorLab_ChromaEstimatedMedian"]
    df["a"] = df["ColorLab_a*Median"]
    df["b"] = df["ColorLab_b*Median"]
    df["L"] = df["ColorLab_L*Median"]
    df["aCV"] = df["ColorLab_a*CoeffVar"]
    df["hue"] = np.degrees(np.arctan2(df.b, df.a)) % 360.0
    df.loc[df.chroma < 1e-6, "hue"] = np.nan
    return df


def main() -> None:
    print("[idea02] loading extract ...")
    df = prep(C.read_extract())
    late = df[(df.tp_h >= LATE[0]) & (df.tp_h <= LATE[1])].copy()
    print(f"  [idea02] late-window rows={len(late)} (passed use non-control plates)")

    # ---------- 2a: species-level pigment triple ----------
    print("[idea02] species x Cu colorimetry ...")
    rows = []
    for (sp, cu), g in late.groupby(["species", "copper_mm"]):
        if pd.isna(sp) or len(g) < 30:
            continue
        cs = circular(g.hue.dropna().values, g.loc[g.hue.notna(), "chroma"].values)
        rows.append({"species": sp, "copper_mm": cu, "n_colonies": len(g),
                     "chroma_median": float(np.nanmedian(g.chroma)),
                     "hue_weighted_mean_deg": cs["mean_deg"],
                     "hue_concentration_R": cs["concentration_R"],
                     "L_median": float(np.nanmedian(g.L))})
    tri = pd.DataFrame(rows)
    C.save(tri, "idea02_species_color_triple.csv")

    # heterogeneity vs Cu per species
    het = (late.groupby(["species", "copper_mm"])["aCV"]
               .agg(["median", "count", "mean"])).reset_index()
    het = het.rename(columns={"median": "aCV_median", "count": "n_colonies", "mean": "aCV_mean"})
    C.save(het, "idea02_heterogeneity.csv")

    # morph clustering on strain x Cu weighted-hue (unit-circle coords)
    print("[idea02] pigment morph clustering ...")
    sc = late.dropna(subset=["hue"]).copy()
    sc["w"] = sc.chroma
    agg = sc.groupby(["strain_id", "copper_mm"]).apply(
        lambda d: pd.Series({
            "hue_wmean": np.degrees(np.arctan2(np.sum(d.w * np.sin(np.radians(d.hue))),
                                               np.sum(d.w * np.cos(np.radians(d.hue))))) % 360,
            "chroma_median": np.nanmedian(d.chroma),
            "n": len(d), "species": d.species.dropna().iloc[0] if d.species.notna().any() else np.nan}),
        include_groups=False).reset_index()
    agg = agg[(agg.n >= 10) & agg.hue_wmean.notna()].copy()
    circ = np.column_stack([np.cos(np.radians(agg.hue_wmean)), np.sin(np.radians(agg.hue_wmean))])
    agg["hue_wmean0"] = agg.hue_wmean.copy()
    best_k, best_sil = 2, -1
    E = {}
    for k in range(2, 7):
        km = KMeans(n_clusters=k, n_init=20, random_state=1).fit(circ)
        if len(np.unique(km.labels_)) < 2:
            continue
        try:
            sil = silhouette_score(circ, km.labels_)
        except Exception:
            sil = -1
        E[k] = (km, sil)
        if sil > best_sil:
            best_sil, best_k = sil, k
    km, sil = E[best_k]
    agg["morph"] = km.labels_
    stab = []
    cents = np.degrees(np.arctan2(km.cluster_centers_[:, 1], km.cluster_centers_[:, 0]))
    cents = np.sort(cents % 360.0)
    for _ in range(50):
        idx = RNG.integers(0, len(circ), len(circ))
        kmb = KMeans(n_clusters=best_k, n_init=10, random_state=1).fit(circ[idx])
        cc = np.degrees(np.arctan2(kmb.cluster_centers_[:, 1], kmb.cluster_centers_[:, 0]))
        cc = np.sort(cc % 360.0)
        d = np.abs(cents - cc)
        stab.append(np.mean(np.minimum(d, 360 - d)))
    morph = agg[["strain_id", "copper_mm", "species", "hue_wmean", "chroma_median", "n", "morph"]]
    C.save(morph, "idea02_pigment_morphs.csv")
    print(f"  [idea02] morph k={best_k} silhouette={best_sil:.3f} centroid-stab={np.mean(stab):.2f} deg")

    # ---------- 2b: run calibration + onset ----------
    print("[idea02] pigment onset + run calibration ...")
    mid = df[(df.tp_h >= 40) & (df.tp_h < 80)].copy()
    ref = mid[mid.species == "Rhodotorula mucilaginosa"]
    cal = (ref.groupby("run_number")[["chroma", "L"]]
              .median().rename(columns={"chroma": "ref_chroma", "L": "ref_L"}))
    cal["n"] = ref.groupby("run_number").size()
    C.save(cal.reset_index(), "idea02_run_calibration.csv")
    base_ch = cal.ref_chroma.median()
    base_L = cal.ref_L.median()
    df = df.merge(cal, on="run_number", how="left")
    df["chroma_adj"] = df.chroma - df.ref_chroma + base_ch
    df["L_adj"] = df.L - df.ref_L + base_L

    # per colony: onset of chroma and darkening
    ons = []
    for (sid, cu), g in df.groupby(["strain_id", "copper_mm"]):
        floor_c = np.nanmedian(g.loc[g.tp_h < 15, "chroma_adj"])
        floor_c = max(floor_c * 2.0, 1.2) if np.isfinite(floor_c) else 1.2
        t_onset, t_dark, pig = np.nan, np.nan, 0
        for wid, w in g.groupby(["run_number", "plate_number", "well_position"]):
            w = w.sort_values("tp_h").drop_duplicates("tp_h")
            if len(w) < 4:
                continue
            c = w.chroma_adj.values
            hit = np.arange(len(w))[c > floor_c]
            for i in hit:
                if i + 1 < len(w) and c[i + 1] > floor_c:
                    t_onset = w.tp_h.iloc[i] if np.isnan(t_onset) else min(t_onset, w.tp_h.iloc[i])
                    pig = 1
                    break
            L = w.L_adj.values
            dL = np.abs(np.gradient(L))
            m = dL > np.nanmedian(dL) * 1.5
            if m.any():
                t = w.tp_h.values[np.where(m)[0][0]]
                t_dark = t if np.isnan(t_dark) else min(t_dark, t)
        ons.append({"strain_id": sid, "copper_mm": cu,
                    "t_chroma_onset_h": t_onset, "t_darkening_h": t_dark,
                    "pigmented": pig,
                    "species": g.species.dropna().iloc[0] if g.species.notna().any() else np.nan})
    onset = pd.DataFrame(ons)
    C.save(onset, "idea02_onset_times.csv")

    # ---------- figures ----------
    plt.figure(figsize=(7, 5))
    for _, r in tri[tri.copper_mm == 0].iterrows():
        plt.scatter(r.chroma_median, r.hue_weighted_mean_deg, s=30, alpha=0.8,
                    label=r.species if r.species in ("Rhodotorula mucilaginosa",
                                                     "Rhodotorula paludigena",
                                                     "Rhodotorula sphaerocarpa",
                                                     "Rhodotorula kratochvilovae") else None)
    plt.xlabel("chroma (pigment amount)"); plt.ylabel("CIELAB hue ° (identity)")
    plt.title("Species pigment identity vs amount, Cu=0 (late window)")
    plt.legend(fontsize=7); plt.tight_layout()
    plt.savefig(C.FIGURES / "fig02a_hue_chroma_species.png", dpi=150)
    plt.close()

    plt.figure(figsize=(6.5, 4.5))
    top = het.species.value_counts().index[:5]
    for sp in top:
        d = het[het.species == sp].sort_values("copper_mm")
        plt.plot(d.copper_mm, d.aCV_median, marker="o", ms=4, label=sp.split()[-1])
    plt.xlabel("copper (mM)"); plt.ylabel("median a*CoeffVar (intra-colony heterogeneity)")
    plt.title("Cu stress and within-colony pigment heterogeneity")
    plt.legend(fontsize=8); plt.tight_layout()
    plt.savefig(C.FIGURES / "fig02a_stress_heterogeneity.png", dpi=150)
    plt.close()

    plt.figure(figsize=(6.5, 5))
    top_m = morph.morph.value_counts().index
    for m in top_m:
        d = morph[morph.morph == m]
        cu0 = d[d.copper_mm == 0]
        plt.scatter(cu0.chroma_median, cu0.hue_wmean, s=14, alpha=0.7,
                    label=f"morph {m} (n={len(cu0)})")
    plt.xlabel("chroma"); plt.ylabel("hue °"); plt.title(f"Pigment morphs (k={best_k}), Cu=0, strains")
    plt.legend(fontsize=7); plt.tight_layout()
    plt.savefig(C.FIGURES / "fig02a_morphs.png", dpi=150)
    plt.close()

    plt.figure(figsize=(7, 5))
    for sp in ("Rhodotorula mucilaginosa", "Rhodotorula paludigena"):
        d = onset[onset.species == sp]
        med = d.t_chroma_onset_h.groupby(d.copper_mm).median()
        plt.plot(med.index, med.values, marker="o", label=sp.split()[-1])
    plt.xlabel("copper (mM)"); plt.ylabel("median pigment onset (h)")
    plt.title("Pigment appearance time vs Cu")
    plt.legend(); plt.tight_layout()
    plt.savefig(C.FIGURES / "fig02b_onset_vs_cu.png", dpi=150)
    plt.close()

    plt.figure(figsize=(5.5, 4))
    plt.bar(cal.reset_index().run_number.astype(str),
            cal.reset_index().ref_chroma.values)
    plt.axhline(base_ch, color="k", ls="--", lw=0.8, label="global median")
    plt.xlabel("imager run"); plt.ylabel("ref chroma (R. mucilaginosa, 40-80h)")
    plt.title("Run-level color calibration offsets"); plt.legend()
    plt.tight_layout(); plt.savefig(C.FIGURES / "fig02b_run_calibration.png", dpi=150)
    plt.close()
    print("[idea02] done.")


def circular(x_deg, w):
    x_deg = np.asarray(x_deg, float)
    w = np.asarray(w, float)
    m = np.isfinite(x_deg)
    x, w = x_deg[m], w[m]
    if len(x) == 0 or w.sum() == 0:
        return {"mean_deg": np.nan, "concentration_R": np.nan, "n": 0}
    t = np.radians(x)
    w = w / w.sum()
    sx, sy = np.sum(w * np.cos(t)), np.sum(w * np.sin(t))
    return {"mean_deg": float(np.degrees(np.arctan2(sy, sx)) % 360),
            "concentration_R": float(np.hypot(sx, sy)), "n": int(m.sum())}


def _ari(a, idx, b):
    """Adjusted Rand index via sklearn between full labels a and bootstrap labels b
    aligned through resampled positions idx."""
    from sklearn.metrics import adjusted_rand_score
    return adjusted_rand_score(a, b)

if __name__ == "__main__":
    sys.exit(main())
