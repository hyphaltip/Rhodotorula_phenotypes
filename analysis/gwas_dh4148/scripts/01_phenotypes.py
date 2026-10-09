#!/usr/bin/env python3
"""Strain-level phenotypes for the GWAS and the dose-0 consistency analysis, from heavy_metal_measurement (curated strain_info).

Outputs (analysis/gwas_dh4148/results/):
  well_growth.csv.gz      per well: relative growth rate (max slope of ln area, 1/h), final window a*, ln area, area at the end
  strain_metal_dose.csv   per strain x metal x dose: mean a*, ln area, RGR over replicate wells
  strain_traits.csv       one row per strain x metal: GWAS traits (a*, area, growth, IC50 ...)
  dose0_wells.csv.gz      dose-0 (YPD) wells of every metal, incl. Fe and Zn (partial plates; used for control consistency only)
Run as a SLURM job from the repo root.
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, duckdb
from scipy.optimize import curve_fit

OUT = Path("analysis/gwas_dh4148/results"); OUT.mkdir(parents=True, exist_ok=True)
METALS = ["Chromium", "Copper", "Lead"]            # D-56: Fe and Zn held out of every dose-response trait
A = "ColorLab_a*GeoMedian"
WIN, W_RGR, MIN_IMG, MIN_AREA = 24.0, 5, 6, 300.0   # a* window (h); images per growth-rate window; images per well; px area at the start of a growth window
log = lambda m: print(m, flush=True)

con = duckdb.connect("db/rhodotorula_phenotypes.duckdb", read_only=True)
d = con.execute(f'''select Metal, run_number, plate_position, Grid_RowNum, Grid_ColNum, Concentration as conc, strain_id, capture_datetime, Shape_Area as area, "{A}" as a
                    from heavy_metal_measurement where strain_id is not null and strain_id not like 'Control%' ''').df()
si = con.execute("select strain_id, sample_name, species, ploidy_status, gwas_panel from strain_info where not is_control").df()
log(f"rows (named strains, all five metals): {len(d):,}")
pk = ["Metal", "run_number", "plate_position"]; wk = pk + ["Grid_RowNum", "Grid_ColNum"]
d["h"] = (d.capture_datetime - d.groupby(pk).capture_datetime.transform("min")).dt.total_seconds() / 3600
assert 90 < d.h.max() < 130
d = d.sort_values("area", ascending=False).drop_duplicates(wk + ["capture_datetime"]).sort_values(wk + ["h"])   # largest object per well-image
d["lnA"] = np.log(d.area)
d["plate_id"] = d.Metal + "_" + d.run_number.astype(str) + "_" + d.plate_position.astype(str)
log(f"well-images: {len(d):,}")

def rgr(g):
    """max OLS slope of ln area against hours over windows of W_RGR consecutive images; window must start at >= MIN_AREA px."""
    h = g.h.values; y = g.lnA.values; ar = g.area.values; best = np.nan
    if len(h) < MIN_IMG: return pd.Series(dict(rgr=np.nan, n_img=len(h)))
    for i in range(len(h) - W_RGR + 1):
        if ar[i] < MIN_AREA or h[i + W_RGR - 1] - h[i] < 4: continue
        s = np.polyfit(h[i:i + W_RGR], y[i:i + W_RGR], 1)[0]
        if np.isnan(best) or s > best: best = s
    return pd.Series(dict(rgr=best, n_img=len(h)))

# growth rate per well
gr = d.groupby(wk).apply(rgr, include_groups=False).reset_index()
# late window values per well (same rule as the a* report: absolute window [T-24, T], T = median plate span per metal; Zinc T = 80 h)
T = d.groupby(pk).h.max().groupby("Metal").median(); T["Zinc"] = 80.0
d["T"] = d.Metal.map(T)
win = d[(d.h >= d["T"] - WIN) & (d.h <= d["T"]) & (d.area >= 2000)]
w = win.groupby(wk).agg(strain_id=("strain_id", "first"), conc=("conc", "first"), plate_id=("plate_id", "first"), n_win=("h", "size"),
                        a=("a", "median"), lnA=("lnA", "median")).reset_index()
w = w[w.n_win >= 2].merge(gr, on=wk, how="left")
w["maxA"] = w.set_index(wk).index.map(d.groupby(wk).area.max())
log(f"wells with window data: {len(w):,}; with growth rate: {w.rgr.notna().sum():,}")
w.to_csv(OUT / "well_growth.csv.gz", index=False)

# dose-0 wells (all five metals) for the control-consistency analysis; growth rate does not need the a* window
g0 = gr.merge(d[d.conc == 0].groupby(wk).agg(strain_id=("strain_id", "first"), plate_id=("plate_id", "first"), a_all=("a", "median"), lnA_all=("lnA", "max")).reset_index(), on=wk)
z0 = w[w.conc == 0][wk + ["a", "lnA"]]
g0 = g0.merge(z0, on=wk, how="left")
g0.to_csv(OUT / "dose0_wells.csv.gz", index=False)
log(f"dose-0 wells: {len(g0):,} ({g0.groupby('Metal').size().to_dict()})")

# strain x metal x dose means (3 metals)
w3 = w[w.Metal.isin(METALS)]
smd = w3.groupby(["Metal", "strain_id", "conc"]).agg(n_wells=("a", "size"), a=("a", "mean"), lnA=("lnA", "mean"), rgr=("rgr", "mean")).reset_index()
smd.to_csv(OUT / "strain_metal_dose.csv", index=False)

def ll(x, ic50, h): return 1.0 / (1.0 + (x / ic50) ** h)
rows = []
for (m, s), g in smd.groupby(["Metal", "strain_id"]):
    g = g.sort_values("conc")
    if g.conc.iloc[0] != 0 or len(g) < 4: continue
    dm = g.conc.max(); x = g.conc.values; xs = x / dm
    r = dict(Metal=m, strain_id=s, n_doses=len(g), a0=g.a.iloc[0], lnA0=g.lnA.iloc[0], rgr0=g.rgr.iloc[0])
    top = g.iloc[-1]
    r.update(a_top=top.a, da_top=top.a - g.a.iloc[0], top_dose=top.conc, top_dose_is_max=bool(top.conc == dm))
    r["a_auc"] = np.trapezoid(g.a.values - g.a.iloc[0], xs)                      # area under the a* change curve (normalised dose 0..1)
    ra = np.exp(g.lnA.values - g.lnA.iloc[0])                                     # relative area (late window), 1 at dose 0
    r["relarea_top"] = ra[-1]; r["relarea_auc"] = np.trapezoid(ra, xs); r["relarea_min"] = ra.min()
    if g.rgr.notna().all() and g.rgr.iloc[0] > 0:
        rr = g.rgr.values / g.rgr.iloc[0]; r["relrgr_auc"] = np.trapezoid(rr, xs); r["relrgr_mid"] = float(np.interp(0.5, xs, rr))
    # IC50 on relative area: reported only when the curve crosses 0.5 inside the tested range and the fit converges
    r["ic50"] = np.nan; r["ic50_status"] = "no_fit"
    if ra.min() >= 0.5: r["ic50_status"] = "right_censored (relative area stays >= 0.5)"
    elif ra[0:2].mean() < 0.5: r["ic50_status"] = "left_censored"
    else:
        try:
            p, _ = curve_fit(ll, x[1:], ra[1:], p0=[np.median(x[x > 0]), 2.0], bounds=([x[x > 0].min() / 5, 0.3], [dm * 5, 12]), maxfev=5000)
            if x[x > 0].min() <= p[0] <= dm: r["ic50"] = p[0]; r["ic50_status"] = "ok"
            else: r["ic50_status"] = "outside_range"
        except Exception as e: r["ic50_status"] = "fit_failed"
    rows.append(r)
tr = pd.DataFrame(rows).merge(si, on="strain_id", how="left")
tr.to_csv(OUT / "strain_traits.csv", index=False)
log(f"strain x metal rows: {len(tr)}; panel rows: {int(tr.gwas_panel.sum())}")
log("IC50 status by metal (all strains):\n" + tr.groupby(["Metal", "ic50_status"]).size().to_string())
log("IC50 status by metal (GWAS panel):\n" + tr[tr.gwas_panel].groupby(["Metal", "ic50_status"]).size().to_string())
