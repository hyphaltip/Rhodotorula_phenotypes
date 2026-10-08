#!/usr/bin/env python3
"""Overview plots of the 5-metal arrayed heavy-metal screen (heavy_metal_measurement).

Run inside a SLURM job from the repo root:
  pixi run python analysis/heavy_metal_overview/scripts/make_overview.py
All decisions are documented in analysis/heavy_metal_overview/HEAVY_METAL_OVERVIEW.md.
Reads Parquet (measurements) and the DB strain_info view (strain table; read-only,
snapshot copy in $SCRATCH if the DB is write-locked). Writes only under
analysis/heavy_metal_overview/results/.
"""
from __future__ import annotations
import pathlib, sys, json
import numpy as np
import pandas as pd
import duckdb
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

HERE = pathlib.Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO / "analysis/gwas/scripts"))
import duckdb_inputs as DI  # read-only DB helper

PARQ = REPO / "data/preprocessed/heavy_metal_array"
RES = HERE.parent / "results"
FIG = RES / "figures"
FIG.mkdir(parents=True, exist_ok=True)

METALS = ["Copper", "Iron", "Lead", "Zinc", "Chromium"]
ABBR = {"Copper": "Cu", "Iron": "Fe", "Lead": "Pb", "Zinc": "Zn", "Chromium": "Cr"}
# Okabe-Ito (colourblind safe)
COL = {"Copper": "#D55E00", "Iron": "#E69F00", "Lead": "#0072B2", "Zinc": "#009E73", "Chromium": "#CC79A7"}
T_MATCH = 90.0       # h since run start; matched timepoint (Pb has no capture 69-87 h; Fe runs end 84-96 h)
T_SENS = 66.0        # sensitivity timepoint
MIN_PAIRED = 30      # min strains paired with 0 dose to accept a 'top dose'
T_TOL = 3.0          # max |capture hour - T_MATCH| accepted per plate
EXPECTED = {"Chromium": 168793, "Copper": 155470, "Iron": 66506, "Lead": 113386, "Zinc": 13216}
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.dpi": 110, "savefig.dpi": 200})


def log(msg: str) -> None:
    print(msg, flush=True)


def conc_label(m: str) -> str:
    return "Concentration (unit not documented; 0-1.2 scale, own axis)" if m == "Chromium" \
        else "Concentration (unit not documented; 0-30 scale)"


# ------------------------------------------------------------------ load
def load():
    con0 = duckdb.connect(":memory:")
    q = f"""SELECT Metal, run_number, plate_position, Grid_RowMajorIdx, capture_datetime,
                   Concentration AS conc, strain_id, Shape_Area AS area,
                   "ColorLab_L*GeoMedian" AS L, "ColorLab_a*GeoMedian" AS a, "ColorLab_b*GeoMedian" AS b
            FROM read_parquet('{PARQ}/*.parquet', union_by_name=true)"""
    d = con0.execute(q).df()
    cnt = d.groupby("Metal").size().to_dict()
    assert cnt == EXPECTED, f"row counts differ from DB (D-35): {cnt}"
    log(f"[load] rows={len(d)} per metal={cnt}")
    con = DI.connect_readonly()
    si = con.execute("SELECT strain_id, is_control, sample_name, species, n_colony_observations FROM strain_info").df()
    log(f"[load] strain_info rows={len(si)} controls={int(si.is_control.sum())}")
    assert len(si) == 334 and si.strain_id.is_unique
    return d, si


# ------------------------------------------------------------------ matched time
def add_hours(d: pd.DataFrame) -> pd.DataFrame:
    t0 = d.groupby(["Metal", "run_number"]).capture_datetime.transform("min")
    d["hours"] = (d.capture_datetime - t0).dt.total_seconds() / 3600.0
    assert d.hours.notna().all() and (d.hours >= 0).all()
    return d


def match_time(d: pd.DataFrame, T: float, tol: float = T_TOL) -> pd.DataFrame:
    """Per plate (Metal, run, plate_position): keep rows of the capture nearest to T h,
    only if |hours - T| <= tol."""
    cap = (d.groupby(["Metal", "run_number", "plate_position", "capture_datetime"]).hours.mean()
             .rename("h").reset_index())
    cap["dev"] = (cap.h - T).abs()
    best = cap.loc[cap.groupby(["Metal", "run_number", "plate_position"]).dev.idxmin()]
    plates_all = len(best)
    ok = best[best.dev <= tol]
    log(f"[match T={T:g}h] plates={plates_all}; within +/-{tol:g} h: {len(ok)}; excluded: {plates_all - len(ok)}")
    excl = best[best.dev > tol].groupby("Metal").size().to_dict()
    log(f"   excluded plates per metal: {excl}; max dev kept={ok.dev.max():.2f} h; median dev={ok.dev.median():.2f} h")
    out = d.merge(ok[["Metal", "run_number", "plate_position", "capture_datetime"]],
                  on=["Metal", "run_number", "plate_position", "capture_datetime"], how="inner")
    log(f"   rows kept={len(out)} of {len(d)}")
    return out, excl, ok


def strain_conc_table(m: pd.DataFrame, si: pd.DataFrame) -> pd.DataFrame:
    """Aggregate colonies (wells) -> median per (Metal, strain, conc). Replicates = wells."""
    n0 = len(m)
    m = m[m.strain_id.notna()]
    log(f"   rows with strain_id: {len(m)} (dropped NULL strain_id: {n0 - len(m)})")
    m = m[m.area.notna()]
    # one well can hold >1 segmented object at one capture; collapse to one value per well (median)
    key = ["Metal", "run_number", "plate_position", "Grid_RowMajorIdx", "strain_id", "conc"]
    nobj = m.groupby(key).size()
    log(f"   wells={len(nobj)}; wells with >1 object={int((nobj > 1).sum())}; objects={len(m)}")
    m = m.groupby(key)[["area", "L", "a", "b"]].median().reset_index()
    g = (m.groupby(["Metal", "strain_id", "conc"])
           .agg(area=("area", "median"), L=("L", "median"), a=("a", "median"), b=("b", "median"),
                n_wells=("area", "size")).reset_index())
    g = g.merge(si[["strain_id", "is_control", "species"]], on="strain_id", how="left")
    assert g.is_control.notna().all(), "strain_id not in strain_info"
    return g


def iqr_band(g, col):
    s = g.groupby("conc")[col]
    return s.median(), s.quantile(0.25), s.quantile(0.75), s.size()


# ------------------------------------------------------------------ plots
def fig_coverage(d, si):
    # (a1) strain x metal: n colony observations (all timepoints), log10
    cnt = d.groupby(["strain_id", "Metal"]).size().unstack("Metal").reindex(columns=METALS)
    nulls = d[d.strain_id.isna()].groupby("Metal").size().reindex(METALS).fillna(0)
    cnt = cnt.reindex(si.strain_id)
    tested = cnt.notna().sum(axis=1)
    ctrl = si.set_index("strain_id").is_control
    order = pd.DataFrame({"t": tested, "c": ctrl.reindex(cnt.index).astype(int),
                          "k": pd.to_numeric(cnt.index, errors="coerce")}).sort_values(["c", "t", "k"], ascending=[True, False, True]).index
    cnt = cnt.loc[order]
    assert len(cnt) == 334
    fig, ax = plt.subplots(1, 2, figsize=(8, 8), gridspec_kw={"width_ratios": [3, 1.2]})
    im = ax[0].imshow(np.log10(cnt.fillna(0).replace(0, np.nan).to_numpy()), aspect="auto",
                      cmap="viridis", interpolation="nearest")
    ax[0].set_xticks(range(5)); ax[0].set_xticklabels([ABBR[m] for m in METALS])
    ax[0].set_ylabel(f"Strain (n={len(cnt)}; 13 Control-N rows at bottom)")
    ax[0].set_yticks([])
    ax[0].set_title("Colony observations per strain and metal\n(all timepoints; blank = not tested)")
    cb = fig.colorbar(im, ax=ax[0], fraction=0.04); cb.set_label("log10(colony observations)")
    sums = d.groupby("Metal").agg(rows=("area", "size"), strains=("strain_id", "nunique"), runs=("run_number", "nunique"))
    sums["null_strain_rows"] = nulls
    sums = sums.reindex(METALS)
    ax[1].barh(range(5), sums.strains, color=[COL[m] for m in METALS])
    ax[1].set_yticks(range(5)); ax[1].set_yticklabels([ABBR[m] for m in METALS]); ax[1].invert_yaxis()
    for i, m in enumerate(METALS):
        ax[1].text(sums.strains[m] + 4, i, f"{int(sums.strains[m])} strains, {int(sums.runs[m])} runs", va="center", fontsize=7)
    ax[1].set_xlim(0, 520); ax[1].set_xlabel("Distinct strain_id (non-NULL)")
    ax[1].set_title("Strains and runs\nper metal")
    fig.suptitle(f"Data coverage: {len(d):,} colony rows; {int(nulls.sum()):,} rows with NULL strain_id not shown by strain", fontsize=9)
    fig.tight_layout(); fig.savefig(FIG / "a1_coverage_strain_by_metal.png"); plt.close(fig)
    sums.to_csv(RES / "coverage_per_metal.csv")
    cnt.to_csv(RES / "coverage_strain_by_metal_counts.csv")
    log(f"[a1] coverage per metal:\n{sums.to_string()}")
    log(f"[a1] strains tested with N metals (non-control): {tested[~ctrl.reindex(tested.index)].value_counts().sort_index().to_dict()}")

    # (a2) runs per metal: distinct strains per run x concentration plate
    fig, axs = plt.subplots(1, 5, figsize=(14, 3.2))
    for ax_, m in zip(axs, METALS):
        s = d[d.Metal == m].groupby(["run_number", "conc"]).strain_id.nunique().unstack("conc")
        im = ax_.imshow(s.to_numpy(), aspect="auto", cmap="cividis", vmin=0, vmax=100)
        ax_.set_xticks(range(s.shape[1])); ax_.set_xticklabels([f"{c:g}" for c in s.columns], fontsize=7)
        ax_.set_yticks(range(s.shape[0])); ax_.set_yticklabels(s.index, fontsize=7)
        for i in range(s.shape[0]):
            for j in range(s.shape[1]):
                v = s.iat[i, j]
                ax_.text(j, i, "" if pd.isna(v) else int(v), ha="center", va="center", color="w", fontsize=6)
        ax_.set_title(f"{ABBR[m]} ({s.shape[0]} runs)")
        ax_.set_xlabel("Concentration" + (" (own scale)" if m == "Chromium" else ""))
    fig.colorbar(im, ax=axs, fraction=0.015, label="distinct strain_id")
    fig.suptitle("Distinct strains per run and concentration (all timepoints; run with 0 = only NULL/control strain rows)", fontsize=9)
    fig.savefig(FIG / "a2_coverage_runs_by_metal.png", bbox_inches="tight"); plt.close(fig)


def fig_dose(g: pd.DataFrame, T: float):
    gs = g[~g.is_control]
    # b: colony size
    fig, axs = plt.subplots(2, 5, figsize=(14, 5.6))
    for j, m in enumerate(METALS):
        x = gs[gs.Metal == m]
        med, q1, q3, n = iqr_band(x, "area")
        axs[0, j].fill_between(med.index, q1, q3, color=COL[m], alpha=0.25, lw=0)
        axs[0, j].plot(med.index, med, "-o", color=COL[m], ms=3)
        axs[0, j].set_title(f"{ABBR[m]}: n={n.min()}-{n.max()} strains/dose")
        axs[0, j].set_xlabel(conc_label(m).replace("Concentration (unit not documented; ", "Conc. (unit?; ").replace(")", ")"), fontsize=7)
        base = x[x.conc == 0].set_index("strain_id").area
        x = x.assign(rel=x.area / x.strain_id.map(base))
        x = x[np.isfinite(x.rel)]
        med, q1, q3, n2 = iqr_band(x, "rel")
        axs[1, j].fill_between(med.index, q1, q3, color=COL[m], alpha=0.25, lw=0)
        axs[1, j].plot(med.index, med, "-o", color=COL[m], ms=3)
        axs[1, j].axhline(1, color="grey", lw=0.6, ls=":")
        axs[1, j].set_xlabel("Conc. (unit?; own scale)" if m == "Chromium" else "Conc. (unit?)", fontsize=7)
        axs[1, j].set_title(f"n={n2.min()}-{n2.max()} strains with 0-dose base", fontsize=8)
    axs[0, 0].set_ylabel("Colony area (Shape_Area, px)\nmedian and IQR across strains")
    axs[1, 0].set_ylabel("Area relative to same strain at 0")
    fig.suptitle(f"Colony size vs concentration at matched time {T:g} h (+/-{T_TOL:g} h per plate); non-control strains; strain value = median over replicate colonies", fontsize=9)
    fig.tight_layout(); fig.savefig(FIG / "b_dose_response_area.png"); plt.close(fig)

    # c: colour
    fig, axs = plt.subplots(3, 5, figsize=(14, 7))
    for i, (col, lab) in enumerate([("L", "L* (GeoMedian)"), ("a", "a* (GeoMedian)"), ("b", "b* (GeoMedian)")]):
        for j, m in enumerate(METALS):
            x = gs[(gs.Metal == m)].dropna(subset=[col])
            med, q1, q3, n = iqr_band(x, col)
            axs[i, j].fill_between(med.index, q1, q3, color=COL[m], alpha=0.25, lw=0)
            axs[i, j].plot(med.index, med, "-o", color=COL[m], ms=3)
            if i == 0:
                axs[i, j].set_title(f"{ABBR[m]}: n={n.min()}-{n.max()} strains/dose")
            if i == 2:
                axs[i, j].set_xlabel("Conc. (unit?; own scale)" if m == "Chromium" else "Conc. (unit?)", fontsize=7)
        axs[i, 0].set_ylabel(f"CIELAB {lab}\nmedian and IQR across strains")
    fig.suptitle(f"Colony colour (CIELAB, per-colony GeoMedian) vs concentration at matched time {T:g} h; non-control strains", fontsize=9)
    fig.tight_layout(); fig.savefig(FIG / "c_dose_response_colour.png"); plt.close(fig)


def tolerance(g: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for m in METALS:
        x = g[(g.Metal == m)]
        piv = x.pivot_table(index="strain_id", columns="conc", values="area")
        # top dose rule: highest dose for which >= MIN_PAIRED strains also have a 0-dose value.
        # (Zn: the 15-30 doses were run on different strains from 0/5/10, so top = 10.)
        paired = {c: int((piv[0.0].notna() & piv[c].notna()).sum()) for c in piv.columns if c > 0}
        top = max(c for c, n in paired.items() if n >= MIN_PAIRED)
        log(f"[tol] {m}: strains paired with 0-dose per dose: {paired} -> top dose {top:g}")
        base, hi = piv.get(0.0), piv.get(top)
        r = pd.DataFrame({"base": base, "top": hi}).dropna()
        r = r[(r.base > 0) & (r.top > 0)]
        r["Metal"] = m; r["top_conc"] = top
        r["log2_ratio"] = np.log2(r.top / r.base)
        log(f"[tol] {m}: strains with colony at 0 and {top:g}: {len(r)}; strains with 0-dose value {base.notna().sum()}, top-dose value {hi.notna().sum()}")
        rows.append(r.reset_index())
    return pd.concat(rows, ignore_index=True)


def fig_cross(tol: pd.DataFrame, T: float):
    w = tol.pivot(index="strain_id", columns="Metal", values="log2_ratio").reindex(columns=METALS)
    rho = pd.DataFrame(index=METALS, columns=METALS, dtype=float); nn = rho.copy()
    for a in METALS:
        for b in METALS:
            ok = w[a].notna() & w[b].notna()
            nn.loc[a, b] = ok.sum()
            rho.loc[a, b] = stats.spearmanr(w.loc[ok, a], w.loc[ok, b])[0] if ok.sum() >= 5 else np.nan
    rho.to_csv(RES / "tolerance_spearman.csv"); nn.to_csv(RES / "tolerance_pairwise_n.csv")
    log(f"[d] Spearman rho:\n{rho.round(3).to_string()}\n[d] pairwise n:\n{nn.astype(int).to_string()}")
    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    im = ax.imshow(rho.astype(float).to_numpy(), cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_xticks(range(5)); ax.set_xticklabels([ABBR[m] for m in METALS])
    ax.set_yticks(range(5)); ax.set_yticklabels([ABBR[m] for m in METALS])
    for i in range(5):
        for j in range(5):
            ax.text(j, i, f"{rho.iat[i, j]:.2f}\nn={int(nn.iat[i, j])}", ha="center", va="center", fontsize=7)
    fig.colorbar(im, label="Spearman rho")
    ax.set_title(f"Cross-metal correlation of tolerance\n(log2 area at top dose / area at 0, {T:g} h)", fontsize=9)
    fig.tight_layout(); fig.savefig(FIG / "d1_tolerance_correlation.png"); plt.close(fig)

    fig, axs = plt.subplots(5, 5, figsize=(10, 10))
    for i, a in enumerate(METALS):
        for j, b in enumerate(METALS):
            ax = axs[i, j]
            if i == j:
                ax.hist(w[a].dropna(), bins=30, color=COL[a])
            else:
                ok = w[a].notna() & w[b].notna()
                ax.scatter(w.loc[ok, b], w.loc[ok, a], s=5, alpha=0.5, color="0.25", lw=0)
            if i == 4: ax.set_xlabel(ABBR[b] + " log2 ratio", fontsize=8)
            if j == 0: ax.set_ylabel(ABBR[a] + " log2 ratio", fontsize=8)
            ax.tick_params(labelsize=6)
    fig.suptitle(f"Per-strain tolerance, log2(area at top dose / area at 0) at {T:g} h; top dose: Cr 1.2, Zn 10, others 30 (unit not documented)", fontsize=9)
    fig.tight_layout(); fig.savefig(FIG / "d2_tolerance_scatter_matrix.png"); plt.close(fig)


def fig_species(tol: pd.DataFrame, si: pd.DataFrame, T: float):
    t = tol.merge(si[["strain_id", "is_control", "species"]], on="strain_id")
    n_in = len(t)
    t = t[~t.is_control & t.species.notna() & (t.species != "Species Not Found")]
    log(f"[e] tolerance rows {n_in} -> {len(t)} after dropping controls / Species Not Found")
    excl = si[(si.species == "Species Not Found") & ~si.is_control]
    log(f"[e] Species Not Found strains (non-control) in strain_info now: {len(excl)}")
    sp_n = t.drop_duplicates("strain_id").species.value_counts()
    order = sp_n.index.tolist()
    fig, axs = plt.subplots(5, 1, figsize=(11, 15), sharex=True)
    rng = np.random.default_rng(7)
    for ax, m in zip(axs, METALS):
        x = t[t.Metal == m]
        nsp = x.groupby("species").size()
        data = [x.loc[x.species == s, "log2_ratio"].to_numpy() for s in order]
        pos = np.arange(len(order))
        bx = [(p, d_) for p, d_ in zip(pos, data) if len(d_) >= 3]
        if bx:
            ax.boxplot([d_ for _, d_ in bx], positions=[p for p, _ in bx], widths=0.6, showfliers=False,
                       medianprops=dict(color="k"), boxprops=dict(color=COL[m]), whiskerprops=dict(color=COL[m]), capprops=dict(color=COL[m]))
        for p, d_ in zip(pos, data):
            ax.scatter(p + rng.uniform(-0.2, 0.2, len(d_)), d_, s=7, color=COL[m], alpha=0.6, lw=0)
        ax.axhline(0, color="grey", lw=0.6, ls=":")
        ax.set_ylabel(f"{ABBR[m]} log2 ratio\n(n={len(x)} strains)")
    axs[-1].set_xticks(np.arange(len(order)))
    axs[-1].set_xticklabels([f"{s.replace('Rhodotorula ', 'R. ')}\n(n={sp_n[s]})" for s in order], rotation=60, ha="right", fontsize=7)
    fig.suptitle(f"Tolerance by species at {T:g} h: log2(area top dose / 0). Boxes for species with >=3 strains with data for that metal; points = strains.\n"
                 f"Excluded: Control-N strains and 'Species Not Found'. n under species = strains in tolerance table (any metal).", fontsize=8)
    fig.tight_layout(); fig.savefig(FIG / "e_tolerance_by_species.png"); plt.close(fig)
    t.groupby(["Metal", "species"]).log2_ratio.agg(["size", "median"]).to_csv(RES / "tolerance_by_species_summary.csv")


def main() -> None:
    d, si = load()
    fig_coverage(d, si)
    d = add_hours(d)
    for m in METALS:
        h = d[d.Metal == m].groupby("run_number").hours.max()
        log(f"[hours] {m}: max hours per run = {h.round(1).to_dict()}")
    sens = {}
    for T in (T_SENS, T_MATCH):
        mt, excl, ok = match_time(d, T)
        g = strain_conc_table(mt, si)
        log(f"   strain x conc cells: {len(g)}")
        sens[T] = (g, tolerance(g))
        if T == T_MATCH:
            ok.to_csv(RES / f"matched_plate_captures_T{int(T)}.csv", index=False)
            g.to_csv(RES / f"strain_conc_medians_T{int(T)}.csv", index=False)
            sens[T][1].to_csv(RES / f"strain_tolerance_T{int(T)}.csv", index=False)
    g, tol = sens[T_MATCH]
    fig_dose(g, T_MATCH)
    fig_cross(tol, T_MATCH)
    fig_species(tol, si, T_MATCH)
    # sensitivity of tolerance to T
    rows = []
    for m in METALS:
        a = sens[T_SENS][1].query("Metal==@m").set_index("strain_id").log2_ratio
        b = sens[T_MATCH][1].query("Metal==@m").set_index("strain_id").log2_ratio
        j = pd.concat([a, b], axis=1, keys=["tS", "tM"]).dropna()
        rows.append({"Metal": m, "n": len(j), "spearman_sens_vs_main": stats.spearmanr(j.tS, j.tM)[0] if len(j) > 4 else np.nan,
                     "median_log2_sens": a.median(), "median_log2_main": b.median()})
    s = pd.DataFrame(rows); s.to_csv(RES / "tolerance_timepoint_sensitivity.csv", index=False)
    log(f"[sens] tolerance T=66 vs 90:\n{s.round(3).to_string(index=False)}")
    log("done")


if __name__ == "__main__":
    main()
