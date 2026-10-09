#!/usr/bin/env python3
"""(1) Which traits correlate with a*?  (2) Do species differ in baseline a* and in stress response?
Same well definition and window as prepare_wells.py (largest object per well-image, window [T-24, T] per metal).
Species come from strain_info (includes user-confirmed overrides); strains labelled 'Species Not Found' are excluded from species tests.
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, duckdb
from scipy import stats
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = Path("analysis/carotenoid_stress_vs_baseline/results"); (OUT / "figures").mkdir(parents=True, exist_ok=True)
WINDOW_H, MIN_IMG, MIN_STRAINS = 24.0, 2, 5
METALS = ["Chromium", "Copper", "Iron", "Lead", "Zinc"]
log = lambda m: print(m, flush=True)
con = duckdb.connect("db/rhodotorula_phenotypes.duckdb", read_only=True)
cols = [r[0] for r in con.execute("describe heavy_metal_measurement").fetchall()
        if r[1] in ("DOUBLE", "BIGINT") and r[0].startswith(("Shape_", "Intensity_", "Texture_", "ColorLab_", "ColorHSV_"))]
log(f"numeric trait columns: {len(cols)}")
sel = ", ".join(f'"{c}"' for c in cols)
df = con.execute(f'''select Metal, run_number, plate_position, Grid_RowNum, Grid_ColNum, Concentration, strain_id, capture_datetime, {sel}
                     from heavy_metal_measurement where strain_id is not null and strain_id not like 'Control%' ''').df()
si = con.execute("select strain_id, species from strain_info").df().set_index("strain_id").species
log(f"rows (named strains): {len(df):,}")
pk = ["Metal", "run_number", "plate_position"]; wk = pk + ["Grid_RowNum", "Grid_ColNum"]
df["h"] = (df.capture_datetime - df.groupby(pk).capture_datetime.transform("min")).dt.total_seconds() / 3600
assert 90 < df.h.max() < 130
df = df.sort_values("Shape_Area", ascending=False).drop_duplicates(wk + ["capture_datetime"])
T = df.groupby(pk).h.max().groupby("Metal").median(); T["Zinc"] = 80.0; df["T"] = df.Metal.map(T)   # Zinc plates stop at 83.7-89.7 h
df = df[df.Shape_Area >= 2000]   # main dataset: largest object >= 2000 px
win = df[(df.h >= df["T"] - WINDOW_H) & (df.h <= df["T"])]
log(f"window rows: {len(win):,}")
agg = {c: "median" for c in cols}; agg.update(strain_id=("strain_id", "first"), conc=("Concentration", "first"), n_img=("capture_datetime", "nunique"))
w = win.groupby(wk).agg(**{c: (c, "median") for c in cols}, strain_id=("strain_id", "first"), conc=("Concentration", "first"), n_img=("capture_datetime", "nunique")).reset_index()
w = w[w.n_img >= MIN_IMG].copy(); w["lnA"] = np.log(w["Shape_Area"]); w["a"] = w["ColorLab_a*GeoMedian"]
w["species"] = w.strain_id.map(si)
log(f"wells: {len(w):,}")
# ---------------- (1) correlations with a*
traits = [c for c in cols if c != "ColorLab_a*GeoMedian" and w[c].notna().mean() > 0.9 and w[c].nunique() > 5]
rows = []
for m in METALS:
    g = w[w.Metal == m]
    for scope, gg in (("all_doses", g), ("dose0", g[g.conc == 0])):
        for c in traits:
            x = gg[[c, "a"]].dropna()
            if len(x) < 50: continue
            rows.append(dict(Metal=m, scope=scope, trait=c, n=len(x), rho=stats.spearmanr(x[c], x.a)[0]))
    # replicate-level: remove strain x dose mean, so only within-strain-within-dose covariation remains
    r = g.copy(); key = ["strain_id", "conc"]
    for c in traits + ["a"]: r[c] = r[c] - r.groupby(key)[c].transform("mean")
    r = r[r.groupby(key).strain_id.transform("size") >= 2]
    for c in traits:
        x = r[[c, "a"]].dropna()
        if len(x) < 50: continue
        rows.append(dict(Metal=m, scope="within_strain_dose", trait=c, n=len(x), rho=stats.spearmanr(x[c], x.a)[0]))
cor = pd.DataFrame(rows); cor.to_csv(OUT / "astar_trait_correlations.csv", index=False)
pd.set_option("display.width", 200); pd.set_option("display.float_format", "{:.2f}".format)
for sc in ("all_doses", "dose0", "within_strain_dose"):
    t = cor[cor.scope == sc].pivot(index="trait", columns="Metal", values="rho")
    t["mean_abs"] = t.abs().mean(axis=1)
    log(f"\n== a* correlation, scope={sc}: top 14 by mean |rho| across metals"); log(t.sort_values("mean_abs", ascending=False).head(14).to_string())
top = cor[cor.scope == "all_doses"].assign(ar=lambda d: d.rho.abs()).groupby("trait").ar.mean().sort_values(ascending=False).head(18).index
piv = cor[(cor.scope == "all_doses") & cor.trait.isin(top)].pivot(index="trait", columns="Metal", values="rho").reindex(top)
fig, ax = plt.subplots(figsize=(7.5, 7)); im = ax.imshow(piv[METALS].values, cmap="RdBu_r", vmin=-1, vmax=1)
ax.set_xticks(range(5)); ax.set_xticklabels(METALS, rotation=45, ha="right"); ax.set_yticks(range(len(piv))); ax.set_yticklabels(piv.index, fontsize=7)
for i in range(piv.shape[0]):
    for j in range(5):
        v = piv[METALS].values[i, j]
        if not np.isnan(v): ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=6)
fig.colorbar(im, label="Spearman rho with a*"); ax.set_title("Traits most correlated with a* (well level, all doses)", fontsize=9)
fig.tight_layout(); fig.savefig(OUT / "figures/fig5_astar_trait_correlations.png", dpi=170); plt.close(fig)
# ---------------- (2) species
ws = w[w.species.notna() & (w.species != "Species Not Found")].copy()
log(f"\nwells with a species label: {len(ws):,} of {len(w):,}")
res = []; box = {}
for m in METALS:
    g = ws[ws.Metal == m]; b = g[g.conc == 0]
    # size-adjusted baseline a*: residual of a ~ lnA on dose-0 wells of this metal
    sl, ic, *_ = stats.linregress(b.lnA, b.a); b = b.assign(a_adj=b.a - (ic + sl * b.lnA))
    sb = b.groupby("strain_id").agg(a=("a", "mean"), a_adj=("a_adj", "mean"), species=("species", "first"), nw=("a", "size"))
    cnt = sb.species.value_counts(); keep = cnt[cnt >= MIN_STRAINS].index
    s = sb[sb.species.isin(keep)]
    for meas in ("a", "a_adj"):
        grp = [x[meas].values for _, x in s.groupby("species")]
        if len(grp) < 2: continue
        H, p = stats.kruskal(*grp); k, n = len(grp), len(s)
        gm = s[meas].mean(); ssb = sum(len(v) * (v.mean() - gm) ** 2 for v in grp); sst = ((s[meas] - gm) ** 2).sum()
        res.append(dict(Metal=m, what=f"baseline_{meas}", n_strains=n, n_species=k, kruskal_H=H, p=p, eta2_H=(H - k + 1) / (n - k), r2_between_species=ssb / sst))
    box[m] = s
    # stress response: change in a* from 0 dose to the highest dose sharing >= 30 strains with 0 dose
    n0 = set(b.strain_id); dmax = max(d for d in g.conc.unique() if d > 0 and len(n0 & set(g[g.conc == d].strain_id)) >= 30)
    top = g[g.conc == dmax].groupby("strain_id").a.mean()
    d = (top - sb.a).dropna().to_frame("delta").join(sb.species); cnt = d.species.value_counts(); d = d[d.species.isin(cnt[cnt >= MIN_STRAINS].index)]
    grp = [x.delta.values for _, x in d.groupby("species")]
    if len(grp) >= 2:
        H, p = stats.kruskal(*grp); k, n = len(grp), len(d)
        res.append(dict(Metal=m, what=f"delta_a_dose{dmax:g}_vs_0", n_strains=n, n_species=k, kruskal_H=H, p=p, eta2_H=(H - k + 1) / (n - k), r2_between_species=np.nan))
R = pd.DataFrame(res); R["p_BH"] = np.nan
idx = R.p.sort_values().index; m_ = len(R); adj = (R.p[idx] * m_ / np.arange(1, m_ + 1)).iloc[::-1].cummin().iloc[::-1].clip(upper=1); R.loc[idx, "p_BH"] = adj.values
R.to_csv(OUT / "species_tests.csv", index=False); log("\n== species tests (species with >= 5 strains; Kruskal-Wallis; eta2_H = effect size)"); log(R.to_string(index=False))
sp_n = pd.concat([x.species.value_counts() for x in box.values()], axis=1); log("\nstrains per species by metal (baseline strains):"); log(sp_n.fillna(0).astype(int).to_string())
fig, ax = plt.subplots(1, 5, figsize=(22, 5.5), sharey=False)
for k, m in enumerate(METALS):
    s = box[m]; order = s.groupby("species").a.median().sort_values().index
    ax[k].boxplot([s[s.species == o].a for o in order], tick_labels=[f"{o.replace('Rhodotorula ','R. ')} ({(s.species==o).sum()})" for o in order], vert=False, showfliers=False)
    ax[k].set_title(f"{m}: baseline a* (0 dose), strain means", fontsize=9); ax[k].tick_params(axis="y", labelsize=7); ax[k].set_xlabel("a*")
fig.suptitle("Baseline a* by species (strain level; species with >= 5 strains; 'Species Not Found' excluded)", fontsize=10)
fig.tight_layout(); fig.savefig(OUT / "figures/fig6_baseline_astar_by_species.png", dpi=170); plt.close(fig)
log("done")
