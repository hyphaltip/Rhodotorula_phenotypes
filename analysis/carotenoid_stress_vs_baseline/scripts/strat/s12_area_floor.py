#!/usr/bin/env python3
"""Where does a* stop depending on colony area? Uses every image (not just the late window): colonies are small early at every dose, so unstressed
small colonies (early) can be compared with stressed small colonies (late) at the same area. Largest object per well-image, named strains only."""
import sys
import numpy as np, pandas as pd, duckdb
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, "analysis/carotenoid_stress_vs_baseline/scripts/strat")
from strat_common import *
F = REP / "figures"; T = REP / "tables"; log = lambda m: print(m, flush=True)
con = duckdb.connect("db/rhodotorula_phenotypes.duckdb", read_only=True)
d = con.execute('''select Metal, run_number, plate_position, Grid_RowNum, Grid_ColNum, Concentration as conc, capture_datetime, Shape_Area as area, "ColorLab_a*GeoMedian" as a
                   from heavy_metal_measurement where strain_id is not null and strain_id not like 'Control%' ''').df()
pk = ["Metal", "run_number", "plate_position"]; wk = pk + ["Grid_RowNum", "Grid_ColNum"]
d["h"] = (d.capture_datetime - d.groupby(pk).capture_datetime.transform("min")).dt.total_seconds() / 3600
d = d.sort_values("area", ascending=False).drop_duplicates(wk + ["capture_datetime"]); d["lnA"] = np.log(d.area)
log(f"well-images: {len(d):,}")
edges = np.log(np.array([100, 200, 400, 700, 1000, 1500, 2000, 3000, 4500, 7000, 10000, 15000, 22000, 33000, 60000])); d["bin"] = pd.cut(d.lnA, edges)
rows = []
for m, g in d.groupby("Metal"):
    dm = g.conc.max()
    for scope, gg in (("dose 0", g[g.conc == 0]), (f"top dose ({dm:g})", g[g.conc == dm]), ("dose 0, first 48 h", g[(g.conc == 0) & (g.h < 48)]), (f"top dose, after 72 h", g[(g.conc == dm) & (g.h > 72)])):
        s = gg.groupby("bin", observed=True).a.agg(["size", "median", lambda x: x.quantile(.25), lambda x: x.quantile(.75), "std"]); s.columns = ["n", "median_a", "q25", "q75", "sd"]; s = s.reset_index()
        s["area_lo"] = s.bin.map(lambda b: float(np.exp(b.left))).astype(float); s["area_hi"] = s.bin.map(lambda b: float(np.exp(b.right))).astype(float); s["Metal"] = m; s["scope"] = scope; rows.append(s.drop(columns="bin"))
tab = pd.concat(rows); tab.to_csv(T / "s12_a_by_area_bin.csv", index=False)
# hinge change point per metal and scope (a ~ lnA, free slope below, different slope above)
hp = []
for m, g in d.groupby("Metal"):
    dm = g.conc.max()
    for scope, gg in (("dose 0", g[g.conc == 0]), ("top dose", g[g.conc == dm]), ("all doses", g)):
        gg = gg[(gg.area > 150) & (gg.area < 40000)]; x, y = gg.lnA.values, gg.a.values; best = None
        for c in np.linspace(np.log(400), np.log(8000), 60):
            if (x < c).sum() < 100 or (x >= c).sum() < 100: continue
            X = np.c_[np.ones_like(x), x, np.maximum(x - c, 0)]; b, *_ = np.linalg.lstsq(X, y, rcond=None); r = ((y - X @ b) ** 2).sum()
            if best is None or r < best[0]: best = (r, c, b)
        if best: r, c, b = best; hp.append(dict(Metal=m, scope=scope, n=len(x), change_point_px=float(np.exp(c)), slope_below=b[1], slope_above=b[1] + b[2], frac_below=float((x < c).mean())))
hp = pd.DataFrame(hp); hp.to_csv(T / "s12_area_hinge.csv", index=False); log(hp.round(2).to_string(index=False))
log("median a* in area bins, dose 0 vs top dose (n>=30):")
for m, g in tab.groupby("Metal"):
    p = g[g.n >= 30].pivot(index="area_lo", columns="scope", values="median_a").round(1); log(m + "\n" + p.to_string())
fig, ax = plt.subplots(1, 5, figsize=(22, 4.6), sharey=True)
for k, m in enumerate(METALS):
    g = tab[tab.Metal == m]
    for scope, col, ls in (("dose 0", "#0072B2", "-"), ("top dose", "#D55E00", "-")):
        z = g[g.scope.str.startswith(scope) & (g.n >= 30) & ~g.scope.str.contains("h")]
        ax[k].plot(np.sqrt(z.area_lo * z.area_hi), z.median_a, "-o", ms=3, color=col, label=scope if scope == "dose 0" else "top dose (all times)")
    z = g[(g.scope == "dose 0, first 48 h") & (g.n >= 30)]; ax[k].plot(np.sqrt(z.area_lo * z.area_hi), z.median_a, "--", color="#0072B2", label="dose 0, first 48 h")
    z = g[(g.scope == "top dose, after 72 h") & (g.n >= 30)]; ax[k].plot(np.sqrt(z.area_lo * z.area_hi), z.median_a, "--", color="#D55E00", label="top dose, after 72 h")
    ax[k].set_xscale("log"); ax[k].axvline(2000, color="k", ls=":", lw=.8); ax[k].set_xlabel("colony area (px, bin geometric mean)"); ax[k].set_title(m, fontsize=10); ax[k].legend(fontsize=6)
ax[0].set_ylabel("median a*")
fig.suptitle("a* against colony area in unstressed and top-dose colonies, all images (dotted line = 2,000 px). Same curve at the same area means a size effect rather than a stress effect", fontsize=9)
fig.tight_layout(); fig.savefig(F / "area_floor_by_metal.png", dpi=160); plt.close(fig); log("done")
