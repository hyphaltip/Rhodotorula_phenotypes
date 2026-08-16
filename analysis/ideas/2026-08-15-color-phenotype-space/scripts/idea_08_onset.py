#!/usr/bin/env python3
"""Idea 08 (Pigment-Onset Dynamics + Growth Coupling).

Is pigmentation gated by colony size (area-dependent switch) or by time?
At well level (Cu=0, typed strains):
  - t_capture = first tp where well-median area > 0.15*max_area
  - t_onset   = first tp where well-median chroma > 2.0  (pigment present)
  - spearman(t_capture, t_onset): area-gating hypothesis predicts tight coupling.
  - pigment-areal pace: slope of log(chroma) vs log(area) over [t_onset, t_onset+30h].
  - logistic onset model (fraction of strains pigmented vs time; Cu 0 vs 5 mM).
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


def main() -> None:
    print("[idea08] loading extract ...")
    df = CM.read_extract()
    d = df[df.strain_id.notna() & ~df.is_control.astype(bool) & (df.tp_h <= 115)].copy()
    d["chr"] = d["ColorLab_ChromaEstimatedMedian"].astype(float)
    d["area"] = d["Shape_Area"].astype(float)
    d["cu"] = np.asarray(d.copper_mm, dtype=float)

    well = (d.groupby(["strain_id", "well_position", "cu"])
               .agg(tp=("tp_h", lambda s: list(s)),
                    chr=("chr", lambda s: list(s)),
                    area=("area", lambda s: list(s)),
                    species=("species", "first")).reset_index())

    rows = []
    for _, r in well.iterrows():
        tp = np.asarray(r.tp, dtype=float)
        ch = np.asarray(r.chr, dtype=float)
        ar = np.asarray(r.area, dtype=float)
        if len(tp) < 6 or ar.max() <= 0:
            continue
        # smooth chronologically
        o = np.argsort(tp); tp, ch, ar = tp[o], ch[o], ar[o]
        med = lambda v: np.array([np.median(v[mask]) for mask in
                                  [np.abs(tp - t) < 3.5 for t in tp]])
        chs, ars = med(ch), med(ar)
        t_cap = next((t for t, a in zip(tp, ars) if a > 0.15 * ars.max()), np.nan)
        t_on = next((t for t, c in zip(tp, chs) if c > 7.0), np.nan)
        rows.append({"strain_id": r.strain_id, "well": r.well_position, "cu": r.cu,
                     "species": r.species, "t_capture_h": t_cap, "t_onset_h": t_on,
                     "max_area": float(ars.max()), "max_chroma": float(chs.max()),
                     "n_tp": int(len(tp))})
    w = pd.DataFrame(rows)

    # --- threshold sensitivity of the onset-capture coupling ---
    sweep = []
    for thr in [2.0, 3.0, 5.0, 7.0, 10.0]:
        on, cap, cu_list = [], [], []
        for _, r in well.iterrows():
            tp = np.asarray(r.tp, dtype=float); ch = np.asarray(r.chr, dtype=float)
            ar = np.asarray(r.area, dtype=float)
            if len(tp) < 6 or ar.max() <= 0:
                continue
            o = np.argsort(tp); tp, ch, ar = tp[o], ch[o], ar[o]
            med = lambda v: np.array([np.median(v[np.abs(tp - t) < 3.5]) for t in tp])
            chs, ars = med(ch), med(ar)
            t_c = next((t for t, a in zip(tp, ars) if a > 0.15 * ars.max()), np.nan)
            t_o = next((t for t, c in zip(tp, chs) if c > thr), np.nan)
            cap.append(t_c); on.append(t_o); cu_list.append(r.cu)
        sw = pd.DataFrame({"t_cap": cap, "t_on": on, "cu": cu_list})
        a = sw[sw.cu == 0]
        pig = a[a.t_on.notna()].dropna(subset=["t_cap"])
        rho = pv = np.nan
        if len(pig) > 50:
            rho, pv = stats.spearmanr(pig.t_cap, pig.t_on)
        sweep.append({"chroma_threshold": thr, "frac_ever_pigment": float(a.t_on.notna().mean()),
                      "median_t_onset_h": float(np.nanmedian(pig.t_on)),
                      "spearman_capture_onset": float(rho), "p": float(pv), "n_pig": len(pig)})
    CM.save(pd.DataFrame(sweep), "idea08_onset_threshold_sweep.csv")
    print("  [idea08] threshold sweep (Cu=0):")
    print(pd.DataFrame(sweep).to_string(index=False))

    # --- A. area-gating: onset vs capture at Cu=0 (wells that pigment) ---
    a0 = w[(w.cu == 0)].copy()
    pig = a0[a0.t_onset_h.notna()].dropna(subset=["t_capture_h"])
    if len(pig) > 20:
        rho, p = stats.spearmanr(pig.t_capture_h, pig.t_onset_h)
        print(f"  [idea08] cu0 pigmenting wells n={len(pig)} "
              f"spearman(t_capture,t_onset)={rho:.3f} p={p:.3g}")
    else:
        rho, p = np.nan, np.nan

    # --- B. pigment-areal pace exponent (log-log coupling after onset) ---
    wt = well.merge(w[["strain_id", "well", "t_onset_h"]], left_on=["strain_id", "well_position"],
                    right_on=["strain_id", "well"], how="left")
    pace_rows = []
    for _, r in wt.iterrows():
        if r.cu != 0 or np.isnan(r.t_onset_h):
            continue
        tp = np.asarray(r.tp, dtype=float); ch = np.asarray(r.chr, dtype=float)
        ar = np.asarray(r.area, dtype=float)
        m = (tp >= r.t_onset_h) & (tp <= r.t_onset_h + 30)
        if m.sum() < 6:
            continue
        lc, la = np.log(ch[m] + 1e-3), np.log(ar[m] + 1)
        if np.std(la) == 0 or np.std(lc) == 0:
            continue
        sl, _, _, _, _ = stats.linregress(la, lc)
        pace_rows.append({"strain_id": r.strain_id, "species": r.species,
                          "pace_loglog": float(sl)})
    pace = pd.DataFrame(pace_rows)
    if len(pace) > 10:
        CM.save(pace, "idea08_pigment_pace.csv")
        print("  [idea08] log-log chroma~area pace by species:")
        print(pace.groupby("species")["pace_loglog"].agg(["count", "mean", "std"])
              .sort_values("mean").to_string())

    # --- C. logistic onset: fraction of strains pigmented vs time ---
    logis = []
    for cu in [0.0, 5.0]:
        sub = w[(w.cu == cu)].copy()
        # one row per strain (worst/earliest per strain: any pigmented colony)
        st = sub.groupby("strain_id").agg(
            t_onset=("t_onset_h", "min"), max_chroma=("max_chroma", "max"),
            species=("species", "first")).reset_index()
        if len(st) < 10:
            continue
        tgrid = np.arange(14, 115, 2.0)
        frac = []
        for t in tgrid:
            frac.append((st.t_onset.notna() & (st.t_onset <= t)).mean())
        # logistic fit: P = 1/(1+exp(-k(t-t50)))
        f = np.asarray(frac)
        if f.max() < 0.15:
            continue
        # linearize t50 at fraction 0.5 (interp)
        from scipy.interpolate import interp1d
        fi = interp1d(f, tgrid, bounds_error=False, fill_value="extrapolate")
        t50 = float(fi(0.5)) if f.min() < 0.5 < f.max() else np.nan
        t25 = float(fi(0.25)) if f.min() < 0.25 < f.max() else np.nan
        logis.append({"cu_mM": cu, "n_strains": len(st), "frac_ever_pigment": float(f.max()),
                      "t25_pigmented_h": t25, "t50_pigmented_h": t50})
        print(f"  [idea08] logistic cu={cu}: frac_ever={f.max():.2f} "
              f"t25={t25 if np.isfinite(t25) else 'NA'} t50={t50 if np.isfinite(t50) else 'NA'}")
    CM.save(pd.DataFrame(logis), "idea08_logistic_onset.csv")

    CM.save(w, "idea08_well_growth_onset.csv")

    # --- figures ---
    plt.figure(figsize=(6.5, 5))
    plt.scatter(pig.t_capture_h, pig.t_onset_h, s=10, alpha=0.5, c="#1f77b4")
    plt.xlabel("time to 15% max colony area (h)")
    plt.ylabel("pigment onset time (h)")
    plt.title(f"Size-gated pigmentation? Cu=0, rho={rho:.2f}")
    if np.isfinite(rho):
        plt.plot([0, 100], [0, 100], "k--", lw=0.8)
    plt.tight_layout(); plt.savefig(CM.FIGURES / "fig08_onset_vs_capture.png", dpi=150)
    plt.close()

    # example curves
    ex = well[(well.cu == 0)].head(2)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    for ax, (_, r) in zip(axes, ex.iterrows()):
        tp = np.asarray(r.tp, dtype=float); o = np.argsort(tp)
        ch = np.polyval(np.polyfit(tp[o], np.asarray(r.chr, dtype=float)[o], 3), tp[o])
        ar = np.asarray(r.area, dtype=float)[o]
        ax.plot(tp[o], ar / ar.max(), label="area / max")
        ax.plot(tp[o], np.clip(ch, 0, None), label="chroma")
        ax.set_xlabel("hours"); ax.legend(fontsize=7); ax.set_title(r.strain_id)
    plt.tight_layout(); plt.savefig(CM.FIGURES / "fig08_example_growth_pigment.png", dpi=150)
    plt.close()
    print("[idea08] done.")


if __name__ == "__main__":
    sys.exit(main())
