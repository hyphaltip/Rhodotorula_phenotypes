#!/usr/bin/env python3
"""Growth curves, relative growth rate, IC50 and trait correlations (figures and tables for the report)."""
from pathlib import Path
import numpy as np, pandas as pd, duckdb
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy import stats
R = Path("analysis/gwas_dh4148/results"); T = Path("analysis/gwas_dh4148/report/tables"); F = Path("analysis/gwas_dh4148/report/figures")
M3 = ["Chromium", "Copper", "Lead"]; COL = {"Chromium": "#0072B2", "Copper": "#D55E00", "Lead": "#CC79A7"}
log = lambda m: print(m, flush=True)
tr = pd.read_csv(R / "strain_traits.csv"); smd = pd.read_csv(R / "strain_metal_dose.csv")
st = pd.read_csv("data/metadata/strain-curation/strain_curation.csv", dtype=str)[["strain_id", "species", "ploidy_status", "gwas_panel", "clade_marker"]].fillna(""); st["strain_id"] = st.strain_id.astype(int)
grp = lambda r: "R. aff. mucilaginosa" if r.clade_marker == "aff_mucilaginosa" else "R. mucilaginosa (haploid)" if (r.species == "Rhodotorula mucilaginosa" and r.ploidy_status == "haploid" and r.clade_marker == "") else ("R. mucilaginosa hybrid diploid" if r.ploidy_status == "diploid_hybrid" and r.species == "Rhodotorula mucilaginosa" else r.species.replace("Rhodotorula ", "R. "))
st["group"] = st.apply(grp, axis=1); tr = tr.drop(columns=["species"]).merge(st[["strain_id", "group"]], on="strain_id")
# 1. growth curves (panel strains, median ln area by hour bin and dose)
con = duckdb.connect("db/rhodotorula_phenotypes.duckdb", read_only=True)
panel = st[st.gwas_panel == "True"].strain_id.tolist()
d = con.execute("select Metal, run_number, plate_position, Grid_RowNum, Grid_ColNum, Concentration as conc, strain_id, capture_datetime, Shape_Area as area from heavy_metal_measurement where Metal in ('Chromium','Copper','Lead') and strain_id is not null and strain_id not like 'Control%'").df()
d["strain_id"] = d.strain_id.astype(int); d = d[d.strain_id.isin(panel)]
pk = ["Metal", "run_number", "plate_position"]; wk = pk + ["Grid_RowNum", "Grid_ColNum"]
d["h"] = (d.capture_datetime - d.groupby(pk).capture_datetime.transform("min")).dt.total_seconds() / 3600
d = d.sort_values("area", ascending=False).drop_duplicates(wk + ["capture_datetime"]); d["lnA"] = np.log(d.area); d["hb"] = (d.h // 6) * 6 + 3
fig, ax = plt.subplots(1, 3, figsize=(15, 4.4), sharey=True)
for k, m in enumerate(M3):
    g = d[d.Metal == m]; doses = sorted(g.conc.unique())
    for i, c in enumerate(doses):
        s = g[g.conc == c].groupby("hb").lnA.agg(["median", "size"]); s = s[s["size"] >= 30]; ax[k].plot(s.index, s["median"], "-o", ms=3, color=plt.cm.viridis(i / max(len(doses) - 1, 1)), label=f"{c:g}")
    ax[k].set_title(f"{m}: median ln area of panel strains", fontsize=9); ax[k].set_xlabel("hours since plate start"); ax[k].legend(fontsize=7, title="dose")
ax[0].set_ylabel("ln colony area (px)"); fig.tight_layout(); fig.savefig(F / "growth_curves_panel.png", dpi=160); plt.close(fig)
# 2. relative growth rate by dose (strain means), panel strains
s = smd[smd.strain_id.isin(panel)].merge(smd[smd.conc == 0][["Metal", "strain_id", "rgr"]].rename(columns={"rgr": "rgr0"}), on=["Metal", "strain_id"])
s["rel"] = s.rgr / s.rgr0; s = s[(s.rgr0 > 0) & s.rel.notna()]
fig, ax = plt.subplots(1, 3, figsize=(15, 4.2), sharey=True); rows = []
for k, m in enumerate(M3):
    g = s[s.Metal == m]; doses = sorted(g.conc.unique()); data = [g[g.conc == c].rel.values for c in doses]
    ax[k].boxplot(data, tick_labels=[f"{c:g}" for c in doses], showfliers=False); ax[k].axhline(1, color="k", lw=.6); ax[k].set_title(f"{m}: growth rate relative to the same strain at dose 0", fontsize=9); ax[k].set_xlabel("dose")
    for c, v in zip(doses, data): rows.append(dict(Metal=m, dose=c, n_strains=len(v), median_rel_rgr=np.median(v), q25=np.quantile(v, .25), q75=np.quantile(v, .75)))
ax[0].set_ylabel("relative maximum growth rate"); fig.tight_layout(); fig.savefig(F / "relative_growth_rate.png", dpi=160); plt.close(fig)
pd.DataFrame(rows).to_csv(T / "relative_growth_rate_by_dose.csv", index=False)
# 3. IC50
ic = tr[tr.ic50_status.isin(["ok"])].copy()
fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
for k, m in enumerate(M3):
    g = tr[tr.Metal == m]; ok = g[g.ic50_status == "ok"]; ax[k].hist(ok.ic50, bins=25, color=COL[m]); cen = (g.ic50_status.str.startswith("right")).sum()
    ax[k].set_title(f"{m}: IC50 of relative area ({len(ok)} fitted; {cen} right-censored; {len(g) - len(ok) - cen} other)", fontsize=8); ax[k].set_xlabel("IC50 (dose units of the screen)")
fig.tight_layout(); fig.savefig(F / "ic50_distribution.png", dpi=160); plt.close(fig)
summ = tr.groupby(["Metal", "ic50_status"]).size().unstack(fill_value=0); summ.to_csv(T / "ic50_status_by_metal.csv")
gs = tr[tr.ic50_status == "ok"].groupby(["Metal", "group"]).ic50.agg(["size", "median", lambda x: x.quantile(.25), lambda x: x.quantile(.75)]); gs.columns = ["n", "median_ic50", "q25", "q75"]; gs = gs[gs.n >= 5].reset_index(); gs.to_csv(T / "ic50_by_group.csv", index=False)
fig, ax = plt.subplots(1, 2, figsize=(13, 5))
for k, m in enumerate(["Chromium", "Lead"]):
    g = tr[(tr.Metal == m) & (tr.ic50_status == "ok")]; gl = [x for x, v in g.groupby("group") if len(v) >= 5]; order = sorted(gl, key=lambda x: g[g.group == x].ic50.median())
    ax[k].boxplot([g[g.group == x].ic50 for x in order], tick_labels=[f"{x} ({(g.group == x).sum()})" for x in order], orientation="horizontal", showfliers=False); ax[k].set_title(f"{m} IC50 by group", fontsize=9); ax[k].tick_params(axis="y", labelsize=7)
fig.tight_layout(); fig.savefig(F / "ic50_by_group.png", dpi=160); plt.close(fig)
# 4. correlations among GWAS traits (panel, raw values)
P = pd.read_csv(R / "gwas/pheno_raw.csv", index_col=0); C = P.corr(method="spearman"); C.round(3).to_csv(T / "gwas_trait_correlations_spearman.csv")
fig, ax = plt.subplots(figsize=(9, 8)); im = ax.imshow(C.values, cmap="RdBu_r", vmin=-1, vmax=1); ax.set_xticks(range(len(C))); ax.set_xticklabels(C.columns, rotation=75, fontsize=7); ax.set_yticks(range(len(C))); ax.set_yticklabels(C.columns, fontsize=7)
fig.colorbar(im, label="Spearman rho (126-strain panel)"); fig.tight_layout(); fig.savefig(F / "gwas_trait_correlations.png", dpi=160); plt.close(fig)
tab = P.describe().T[["count", "mean", "std", "min", "50%", "max"]].reset_index().rename(columns={"index": "trait"}); tab.to_csv(T / "gwas_trait_summary.csv", index=False)
log(summ.to_string()); log(gs.round(3).to_string(index=False))
