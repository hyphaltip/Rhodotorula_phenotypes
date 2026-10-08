#!/usr/bin/env python3
"""Zinc rescue. Zinc has 9 plates (one per dose) in two runs with disjoint strain sets: set A (run d000388: doses 0,5,10,15) and set B (run d000390: doses 10-30).
Because a plate effect is the same for every strain on that plate, anything based on differences BETWEEN strains within a plate (ranks, deviations from the plate mean,
correlations across strains) is not confounded with plate. The population-mean dose effect is (one plate per dose)."""
import sys, itertools
import numpy as np, pandas as pd, duckdb
from scipy import stats
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, "analysis/carotenoid_stress_vs_baseline/scripts/strat")
from strat_common import *
import os
MINAREA = float(os.environ.get("ZN_MINAREA", "2000")); TAG = "" if MINAREA > 0 else "_nofilter"   # main = largest object >= 2000 px; ZN_MINAREA=0 gives the unfiltered comparison
F = REP / "figures"; T = REP / "tables"; log = lambda m: print(m, flush=True)
con = duckdb.connect("db/rhodotorula_phenotypes.duckdb", read_only=True)
z = con.execute('''select run_number, plate_position, Grid_RowNum, Grid_ColNum, Concentration as conc, strain_id, capture_datetime, Shape_Area as area, "ColorLab_a*GeoMedian" as a, "ColorLab_b*GeoMedian" as b
                   from heavy_metal_measurement where Metal='Zinc' and strain_id is not null and strain_id not like 'Control%' ''').df()
pk = ["run_number", "plate_position"]; wk = pk + ["Grid_RowNum", "Grid_ColNum"]
z["h"] = (z.capture_datetime - z.groupby(pk).capture_datetime.transform("min")).dt.total_seconds() / 3600
z = z.sort_values("area", ascending=False).drop_duplicates(wk + ["capture_datetime"]); z["lnA"] = np.log(z.area)
if MINAREA > 0: z = z[z.area >= MINAREA]   # same area cutoff as the main dataset
log(f"min area {MINAREA:g} px; well-images: {len(z):,}")
span = z.groupby(pk).agg(conc=("conc", "first"), imgs=("capture_datetime", "nunique"), span_h=("h", "max")).reset_index(); log("plates:\n" + span.to_string(index=False))
def wells(Tend, width=24.0):
    w = z[(z.h >= Tend - width) & (z.h <= Tend)].groupby(wk).agg(strain_id=("strain_id", "first"), conc=("conc", "first"), n_img=("capture_datetime", "nunique"), a=("a", "median"), b=("b", "median"), lnA=("lnA", "median")).reset_index()
    return w[w.n_img >= 2].copy()
cover = []
for Tend in (66, 80, 90, 107.6):
    w = wells(Tend); c = w.groupby(pk).size().reindex(pd.MultiIndex.from_frame(span[pk]), fill_value=0).values; cover.append(pd.Series(c, name=f"T={Tend:g}h"))
cov = pd.concat([span.reset_index(drop=True), pd.concat(cover, axis=1)], axis=1); cov.to_csv(T / f"s11_zinc_window_coverage{TAG}.csv", index=False); log("wells per plate by window end:\n" + cov.to_string(index=False))
w = wells(80.0); w["set"] = np.where(w.run_number == "d000388", "A (run 388)", "B (run 390)"); w["plate"] = w.run_number + "_" + w.plate_position.astype(str)
w.to_csv(R / f"zinc_wells_T80{TAG}.csv", index=False); log(f"T=80 h window [56,80]: {len(w)} wells, plates {w.plate.nunique()}, strains {w.strain_id.nunique()}")
st = pd.read_csv(R / "strat/strain_table.csv", dtype={"strain_id": str}); w["strain_id"] = w.strain_id.astype(str); w = w.merge(st[["strain_id", "species"]], on="strain_id", how="left")
log("species per set (strains):\n" + w.groupby("set").apply(lambda g: g.drop_duplicates("strain_id").species.value_counts().head(4).to_dict()).to_string())
# ---- set A paired analysis
A = w[w.set.str.startswith("A")]; dA = sorted(A.conc.unique()); B = w[w.set.str.startswith("B")]; dB = sorted(B.conc.unique())
log(f"set A doses {dA}, strains per dose {A.groupby('conc').strain_id.nunique().to_dict()}; set B doses {dB}, strains per dose {B.groupby('conc').strain_id.nunique().to_dict()}")
rows = []
for nm, g, ds in (("A", A, dA), ("B", B, dB)):
    for col in ("a", "lnA", "b"):
        P = g.pivot_table(index="strain_id", columns="conc", values=col)
        for d1, d2 in itertools.combinations(ds, 2):
            x = P[[d1, d2]].dropna()
            if len(x) >= 10: rows.append(dict(set=nm, trait=col, dose1=d1, dose2=d2, n_strains=len(x), spearman_between_plates=stats.spearmanr(x[d1], x[d2])[0], mean_diff=(x[d2] - x[d1]).mean(), sd_diff=(x[d2] - x[d1]).std()))
rk = pd.DataFrame(rows); rk.to_csv(T / f"s11_zinc_rank_consistency{TAG}.csv", index=False); log(rk[rk.trait == "a"].round(2).to_string(index=False)); log(rk[rk.trait == "lnA"].round(2).to_string(index=False))
# ---- set exchangeability at the doses both sets share
sh = []
for d in sorted(set(dA) & set(dB)):
    for col in ("a", "lnA"):
        xa, xb = A[A.conc == d][col], B[B.conc == d][col]
        if len(xa) >= 10 and len(xb) >= 10: sh.append(dict(dose=d, trait=col, n_A=len(xa), n_B=len(xb), mean_A=xa.mean(), mean_B=xb.mean(), mann_whitney_p=stats.mannwhitneyu(xa, xb).pvalue))
sh = pd.DataFrame(sh); sh.to_csv(T / f"s11_zinc_sets_at_shared_doses{TAG}.csv", index=False); log("sets at shared doses:\n" + sh.round(3).to_string(index=False))
# ---- cross-metal: zinc strain deviations (within plate) vs other metals' strain-level response
o = pd.read_csv(R / "wells.csv"); o["strain_id"] = o.strain_id.astype(str); o = o[o.Metal != "Zinc"]
resp = {}
for m, g in o.groupby("Metal"):
    dm = g.conc.max(); s = g[g.conc.isin([0, dm])].groupby(["strain_id", "conc"])[["a", "lnA"]].mean().unstack("conc").dropna()
    resp[m] = pd.DataFrame({f"{m} d_a": s[("a", dm)] - s[("a", 0.0)], f"{m} d_lnA": s[("lnA", dm)] - s[("lnA", 0.0)], f"{m} a_top": s[("a", dm)], f"{m} a_0": s[("a", 0.0)]})
R_ = pd.concat(resp.values(), axis=1)
cm = []
for nm, g, d0, d1 in (("A", A, 0.0, 15.0), ("B", B, 15.0, 30.0)):
    P = g.pivot_table(index="strain_id", columns="conc", values=["a", "lnA"])
    zs = pd.DataFrame({"Zn a at dose %g" % d1: P[("a", d1)], "Zn d_a (%g vs %g)" % (d1, d0): P[("a", d1)] - P[("a", d0)], "Zn d_lnA (%g vs %g)" % (d1, d0): P[("lnA", d1)] - P[("lnA", d0)], "Zn lnA at dose %g" % d1: P[("lnA", d1)]})
    j = zs.join(R_, how="inner")
    for zc in zs.columns:
        for oc in R_.columns:
            x = j[[zc, oc]].dropna()
            if len(x) >= 20: cm.append(dict(zn_set=nm, zinc_measure=zc, other=oc, n_strains=len(x), spearman=stats.spearmanr(x[zc], x[oc])[0]))
cm = pd.DataFrame(cm); cm.to_csv(T / f"s11_zinc_cross_metal{TAG}.csv", index=False)
log("cross-metal: strongest |rho| per zinc measure"); log(cm.assign(ar=cm.spearman.abs()).sort_values("ar", ascending=False).groupby(["zn_set", "zinc_measure"]).head(3)[["zn_set", "zinc_measure", "other", "n_strains", "spearman"]].round(2).to_string(index=False))
# ---- figures
fig, ax = plt.subplots(1, 2, figsize=(14, 4.6)); x = np.arange(len(cov)); lab = [f"{r.run_number[-3:]}/{r.plate_position}\ndose {r.conc:g}" for r in cov.itertuples()]
for j_, (c, col) in enumerate(zip(["T=107.6h", "T=80h"], ["#999999", "#0072B2"])): ax[0].bar(x + (j_ - .5) * .38, cov[c], .36, color=col, label=f"window ends {c[2:]}")
ax[0].set_xticks(x); ax[0].set_xticklabels(lab, fontsize=7); ax[0].set_ylabel("wells usable (>= 2 images in a 24 h window)"); ax[0].legend(); ax[0].set_title("Zinc: wells per plate by late-window end", fontsize=9)
ax[1].bar(x, cov.span_h, color="#D55E00"); ax[1].axhline(80, color="k", ls="--", lw=.8); ax[1].set_xticks(x); ax[1].set_xticklabels(lab, fontsize=7); ax[1].set_ylabel("imaging span (h)"); ax[1].set_title("Plate imaging span (dashed = 80 h)", fontsize=9)
fig.tight_layout(); fig.savefig(F / f"zn_window_coverage{TAG}.png", dpi=160); plt.close(fig)
fig, ax = plt.subplots(2, 3, figsize=(16, 8))
for r_, (nm, g, ds) in enumerate((("A (run 388, 80 strains)", A, dA), ("B (run 390, about 70 strains)", B, dB))):
    for c_, (col, lab_) in enumerate((("a", "a*"), ("lnA", "ln colony area"), ("b", "b*"))):
        P = g.pivot_table(index="strain_id", columns="conc", values=col)
        for _, row in P.iterrows(): ax[r_, c_].plot(P.columns, row.values, color="#999999", lw=.3, alpha=.4)
        m_ = P.mean(); s_ = P.sem(); ax[r_, c_].errorbar(P.columns, m_, yerr=1.96 * s_, fmt="-o", color="#0072B2", lw=2, capsize=3); ax[r_, c_].set_xlabel("Zn conc"); ax[r_, c_].set_ylabel(lab_); ax[r_, c_].set_title(f"Set {nm}: {lab_} (thin = single strains)", fontsize=8)
fig.suptitle("Zinc, window [56, 80] h: strains followed across doses within a set. The mean line also contains the plate effect (one plate per dose)", fontsize=10); fig.tight_layout(); fig.savefig(F / f"zn_paired_dose_response{TAG}.png", dpi=160); plt.close(fig)
fig, ax = plt.subplots(2, 4, figsize=(18, 8))
for r_, (nm, g, ds) in enumerate((("A", A, dA), ("B", B, dB))):
    for c_, col in enumerate(("a", "lnA")):
        P = g.pivot_table(index="strain_id", columns="conc", values=col); C_ = P.corr(method="spearman"); im = ax[r_, c_].imshow(C_.values, cmap="RdBu_r", vmin=-1, vmax=1)
        ax[r_, c_].set_xticks(range(len(ds))); ax[r_, c_].set_xticklabels([f"{d:g}" for d in P.columns]); ax[r_, c_].set_yticks(range(len(ds))); ax[r_, c_].set_yticklabels([f"{d:g}" for d in P.columns])
        for i in range(C_.shape[0]):
            for j in range(C_.shape[1]): ax[r_, c_].text(j, i, f"{C_.values[i, j]:.2f}", ha="center", va="center", fontsize=8)
        ax[r_, c_].set_title(f"Set {nm}: strain rank agreement between dose plates, {'a*' if col == 'a' else 'ln area'}", fontsize=8)
sub = cm[cm.zinc_measure.str.contains("d_lnA") | cm.zinc_measure.str.contains("d_a")]
for r_, nm in enumerate(("A", "B")):
    for c_, meas in enumerate(("d_lnA", "d_a")):
        s = sub[(sub.zn_set == nm) & sub.zinc_measure.str.contains(meas)]; piv = s[s.other.str.contains(meas)].set_index("other").spearman
        ax[r_, 2 + c_].bar(range(len(piv)), piv.values, color="#009E73"); ax[r_, 2 + c_].set_xticks(range(len(piv))); ax[r_, 2 + c_].set_xticklabels([o.replace(" " + meas, "") for o in piv.index], rotation=30, fontsize=8)
        ax[r_, 2 + c_].axhline(0, color="k", lw=.5); ax[r_, 2 + c_].set_ylim(-.6, .8); ax[r_, 2 + c_].set_title(f"Set {nm}: Zn {meas} vs same strains' {meas} in other metals (Spearman)", fontsize=8)
fig.suptitle("Zinc: within-plate strain comparisons are not affected by the plate effect", fontsize=10); fig.tight_layout(); fig.savefig(F / f"zn_strain_level_consistency{TAG}.png", dpi=160); plt.close(fig)

# ---- composite dose series: set A (0-15) and set B (15-30) bridged at the dose both sets share (15)
comp = []
for nm, g in (("A", A), ("B", B)):
    for d, gg in g.groupby("conc"):
        r_ = dict(set=nm, dose=d, n_strains=gg.strain_id.nunique())
        for col in ("a", "b", "lnA"):
            sm = gg.groupby("strain_id")[col].mean(); r_[f"{col}_mean"] = sm.mean(); r_[f"{col}_ci95"] = 1.96 * sm.sem()
        comp.append(r_)
comp = pd.DataFrame(comp); comp.to_csv(T / f"s11_zinc_composite_series{TAG}.csv", index=False); log("composite series:\n" + comp.round(2).to_string(index=False))
fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
for c_, (col, lab_) in enumerate((("a", "a*"), ("b", "b*"), ("lnA", "ln colony area"))):
    for nm, colr in (("A", "#0072B2"), ("B", "#D55E00")):
        z_ = comp[comp.set == nm]; ax[c_].errorbar(z_.dose, z_[f"{col}_mean"], yerr=z_[f"{col}_ci95"], fmt="-o", color=colr, capsize=3, label=f"set {nm}: run {'388' if nm == 'A' else '390'}, {int(z_.n_strains.min())}-{int(z_.n_strains.max())} strains")
    ax[c_].axvline(15, color="k", ls=":", lw=.8); ax[c_].set_xlabel("Zn conc"); ax[c_].set_ylabel(lab_); ax[c_].legend(fontsize=7); ax[c_].set_title(f"Zinc composite: {lab_} (strain means, 95% CI)", fontsize=9)
fig.suptitle("Zinc as one experiment: set A (doses 0-15) and set B (doses 15-30) agree at the shared dose 15 (dotted line); one plate per dose, so the level at each dose includes a plate effect", fontsize=9)
fig.tight_layout(); fig.savefig(F / f"zn_composite_series{TAG}.png", dpi=160); plt.close(fig)
log("done")
