#!/usr/bin/env python3
"""Assemble report/REPORT.md (tables from CSV, PNGs embedded with relative paths so GitHub renders them) and copy figures/tables."""
import shutil, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, "analysis/carotenoid_stress_vs_baseline/scripts/strat")
from strat_common import *
T = REP / "tables"; F = REP / "figures"
for f in ("fig1_astar_vs_size_by_dose", "fig2_baseline_vs_stressed", "fig3_strain_baseline_vs_induction", "fig4_repeatability", "fig5_astar_trait_correlations", "fig6_baseline_astar_by_species"):
    shutil.copy(R / "figures" / (f + ".png"), F / (f + ".png"))
for f in ("mixed_model_summary", "dose_factor_effects_M2", "astar_trait_correlations", "species_tests"):
    shutil.copy(R / (f + ".csv"), T / (f + ".csv"))
def md(df, floatfmt="{:.2f}"):
    df = df.copy(); df.columns = [str(c) for c in df.columns]; cols = list(df.columns)
    for c in cols:
        if df[c].dtype.kind == "f": df[c] = df[c].map(lambda v: "" if pd.isna(v) else ("<0.001" if (c.lower().startswith("p") and 0 <= v < 0.001) else floatfmt.format(v)))
    esc = lambda x: str(x).replace("|", "\\|").replace("*", "\\*")   # pipes and asterisks break GitHub tables
    rows = [[esc(x) for x in r] for r in df.itertuples(index=False)]
    wd = [min(40, max(len(c), *(len(r[i]) for r in rows))) if rows else len(c) for i, c in enumerate(cols)]  # proportional separators -> pandoc column widths
    out = ["| " + " | ".join(esc(c) for c in cols) + " |", "|" + "|".join("-" * max(3, w) for w in wd) + "|"]
    out += ["| " + " | ".join(r) + " |" for r in rows]
    return "\n".join(out)
img = lambda f, cap: "![](figures/%s.png)\n\n*%s*\n" % (f, cap.replace("\\*", "*").replace("*", "\\*"))   # empty alt text: pandoc would otherwise add a second caption in the PDF

# ---- fig5 redrawn from the saved correlation table (rotated x labels)
_c = pd.read_csv(T / "astar_trait_correlations.csv"); _c = _c[_c.scope == "all_doses"]
_top = _c.assign(ar=_c.rho.abs()).groupby("trait").ar.mean().sort_values(ascending=False).head(18).index
_p = _c[_c.trait.isin(_top)].pivot(index="trait", columns="Metal", values="rho").reindex(_top)[METALS]
fig, ax = plt.subplots(figsize=(8.2, 7.4)); im = ax.imshow(_p.values, cmap="RdBu_r", vmin=-1, vmax=1)
ax.set_xticks(range(5)); ax.set_xticklabels(METALS, rotation=45, ha="right", fontsize=9); ax.set_yticks(range(len(_p))); ax.set_yticklabels(_p.index, fontsize=7)
for i in range(_p.shape[0]):
    for j in range(5):
        v = _p.values[i, j]
        if not np.isnan(v): ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=6)
fig.colorbar(im, label="Spearman rho with a*"); ax.set_title("Traits most correlated with a* (well level, all doses)", fontsize=9)
fig.tight_layout(); fig.savefig(F / "fig5_astar_trait_correlations.png", dpi=170); plt.close(fig)
# phylo figure: lambda only (Blomberg K is unstable here, see text)
ph = pd.read_csv(T / "s2_phylo_signal.csv"); sub = ph[ph.trait.str.startswith(("baseline a*", "change"))].copy()
fig, ax = plt.subplots(figsize=(12, 4.6)); labs = []
for j, sc in enumerate(("all_strains_with_tip", "R_mucilaginosa_only")):
    s = sub[sub.scope == sc]; x = np.arange(len(s)) + (j - .5) * .38; ax.bar(x, s.pagel_lambda, .36, color=("#0072B2", "#D55E00")[j], label=sc.replace("_", " "))
labs = [m[:2] + " " + ("base" if t == "baseline a*" else "base adj" if "size" in t else "change") for m, t in zip(s.Metal, s.trait)]
ax.set_xticks(range(len(s))); ax.set_xticklabels(labs, rotation=70, fontsize=8); ax.set_ylabel("Pagel's lambda"); ax.axhline(0, color="k", lw=.5); ax.legend(fontsize=8)
ax.set_title("Phylogenetic signal (Pagel's lambda) of baseline a* and change in a* at top dose", fontsize=10); fig.tight_layout(); fig.savefig(F / "s2_phylo_signal.png", dpi=170); plt.close(fig)
# ---- tables
mm = pd.read_csv(T / "mixed_model_summary.csv"); mm["model"] = mm.model.map({"M0": "total", "M1": "at fixed size"})
t_mm = md(mm[["Metal", "model", "n_wells", "n_strains", "dose_effect", "dose_se", "dose_p", "repeatability_dose0", "repeatability_dose1", "singular"]].rename(columns={"dose_effect": "dose effect (0 to top)", "dose_se": "SE", "dose_p": "p", "repeatability_dose0": "strain share @0", "repeatability_dose1": "strain share @top"}))
m2 = pd.read_csv(T / "dose_factor_effects_M2.csv"); t_m2 = md(m2.pivot(index="conc", columns="Metal", values="effect_at_fixed_size").reset_index().rename(columns={"conc": "dose"}).fillna(np.nan))
cor = pd.read_csv(T / "astar_trait_correlations.csv"); cor = cor[cor.trait != "ColorLab_a*Medoid"]
def topcor(sc, n=10):
    p = cor[cor.scope == sc].pivot(index="trait", columns="Metal", values="rho"); p["mean abs rho"] = p.abs().mean(axis=1); return md(p.sort_values("mean abs rho", ascending=False).head(n).reset_index())
om = pd.read_csv(T / "s1_omnibus_tests.csv"); t_om = md(om[["Metal", "model", "F", "p", "p_BH"]])
sb = pd.read_csv(T / "s1_species_baseline.csv"); sb = sb[sb.size_adjusted == True]; sb["species"] = sb.species.str.replace("Rhodotorula ", "R. ")
t_sb = md(sb[["Metal", "species", "n_strains", "estimate", "se", "p", "p_BH"]])
ss = pd.read_csv(T / "s1_species_dose_slope.csv"); ss["species"] = ss.species.str.replace("Rhodotorula ", "R. "); t_ss = md(ss[["Metal", "species", "n_strains", "slope_full_range", "se", "vs_reference_p"]].rename(columns={"vs_reference_p": "p vs mucilaginosa"}))
po = pd.read_csv(T / "s1_pooled_other_vs_mucilaginosa.csv"); t_po = md(po)
pb = pd.read_csv(T / "s1_population_baseline.csv"); pb["population"] = pb.population.str.replace("poppop", "pop"); t_pb = md(pb)
pp = pd.read_csv(T / "s1_population_dose_slope.csv"); t_pp = md(pp)
t_ph = md(sub[["Metal", "scope", "trait", "n_strains", "pagel_lambda", "p_lambda_gt_0"]].rename(columns={"p_lambda_gt_0": "p (lambda > 0)"}))
sr = pd.read_csv(T / "s3_size_ratio_by_dose.csv"); t_sr = md(sr[["Metal", "conc", "median_ratio", "q25", "q75", "n_strains", "regime"]])
rs = pd.read_csv(T / "s3_regime_slopes.csv"); t_rs = md(rs)
s4 = pd.read_csv(T / "s4_size_matched.csv"); t_s4 = md(s4[["Metal", "size_bin", "n_wells", "n_doses", "mean_a", "sd_a", "slope_full_range", "se", "p", "max_dose_s_in_bin", "change_over_observed_range"]])
s5 = pd.read_csv(T / "s5_split_half.csv"); s5 = s5[s5.measure == "a_adj"]
t_s5 = md(s5[["Metal", "top_dose", "n_strains", "naive_gap", "naive_ci_lo", "naive_ci_hi", "split_gap_mean", "split_gap_ci_lo", "split_gap_ci_hi", "reliability_change_A_vs_B", "reliability_baseline_A_vs_B"]])
h6 = pd.read_csv(T / "s6_run_heterogeneity.csv"); t_h6 = md(h6); p6 = pd.read_csv(T / "s6_per_run_dose_effect.csv"); t_p6 = md(p6); v6 = pd.read_csv(T / "s6_variance_components.csv"); t_v6 = md(v6[["Metal", "n_runs", "share_strain", "share_run", "share_plate", "share_resid", "singular"]])
sj = pd.read_csv(T / "s2_jitter_sensitivity.csv"); t_sj = md(sj[["Metal", "jitter", "cond", "lam", "K", "pK"]], "{:.3g}")

# ---- new sections: retained wells, area floor, species strata, b*/morphology, Zinc, sensitivity
RS = {}
nf = pd.read_csv(R / "wells_nofilter.csv"); mn = pd.read_csv(R / "wells.csv")
a_ = nf.groupby(["Metal", "conc"]).size().rename("unfiltered"); b_ = mn.groupby(["Metal", "conc"]).size().rename("main"); ret = pd.concat([a_, b_], axis=1).fillna(0).reset_index()
ret["retained_%"] = 100 * ret.main / ret.unfiltered; ret.to_csv(T / "wells_retained_by_dose.csv", index=False)
ret["dose"] = ret.conc.map(lambda x: f"{x:g}"); RS["T_RETAINED"] = md(ret[["Metal", "dose", "unfiltered", "main", "retained_%"]].rename(columns={"unfiltered": "wells, unfiltered", "main": "wells, main"}), "{:.0f}")
tot = ret.groupby("Metal")[["unfiltered", "main"]].sum(); RS["T_RETAINED_NOTE"] = "Share of unfiltered wells kept: " + ", ".join(f"{m} {100 * r.main / r.unfiltered:.0f}%" for m, r in tot.iterrows()) + "."
bins = pd.read_csv(T / "s12_a_by_area_bin.csv"); b0 = bins[(bins.scope == "dose 0") & (bins.n >= 30)].copy(); b0["area (px, bin start)"] = b0.area_lo.round(-1).astype(int)
RS["T_S12_BINS"] = md(b0.pivot(index="area (px, bin start)", columns="Metal", values="median_a").reset_index(), "{:.1f}")
hg = pd.read_csv(T / "s12_area_hinge.csv"); RS["T_S12_HINGE"] = md(hg[["Metal", "scope", "n", "change_point_px", "slope_below", "slope_above", "frac_below"]].rename(columns={"change_point_px": "change-point (px)", "frac_below": "share of wells below"}), "{:.2f}").replace(".00 |", " |")
RS["T_S9_COUNTS"] = md(pd.read_csv(T / "s9_species_strain_counts.csv"))
sm_ = pd.read_csv(T / "s9_stratum_models.csv"); sm_["stratum"] = sm_.stratum.str.replace("Rhodotorula ", "R. ")
pv = sm_[sm_.model == "total"].pivot(index="stratum", columns="Metal", values="dose_effect")
RS["T_S9_MODEL"] = md(pv.reset_index(), "{:.1f}")
s10 = pd.read_csv(T / "s10_trait_dose_effects.csv"); s10 = s10[s10.dataset.str.startswith("main")]; ordr = ["L", "a", "b", "chroma", "hue_deg", "sat", "val", "circ", "solid", "ecc", "compact", "extent", "aspect"]
for k_, mod in (("T_S10_TOTAL", "total"), ("T_S10_FIXED", "at fixed size")):
    z_ = s10[s10.model == mod].pivot(index="trait", columns="Metal", values="effect_sd").reindex(ordr); RS[k_] = md(z_.reset_index(), "{:.1f}")
zc = pd.read_csv(T / "s11_zinc_window_coverage.csv"); zc["plate"] = zc.run_number + "/" + zc.plate_position.astype(str); RS["T_ZN_COVER"] = md(zc[["plate", "conc", "imgs", "span_h", "T=66h", "T=80h", "T=90h", "T=107.6h"]].rename(columns={"conc": "dose", "imgs": "images", "span_h": "span (h)"}), "{:.1f}")
zp = pd.read_csv(T / "s11_zinc_composite_series.csv"); RS["T_ZN_COMP"] = md(zp, "{:.2f}")
zn_ = pd.read_csv(T / "s11_zinc_composite_series_nofilter.csv"); RS["T_ZN_COMP_NF"] = md(zn_, "{:.2f}")
cu_ = pd.read_csv(R / "wells_traits.csv"); cu_ = cu_[cu_.Metal == "Copper"].groupby(["strain_id", "conc"])[["a", "b", "lnA"]].mean().groupby("conc").mean().reset_index(); cu_.to_csv(T / "s10_cu_dose_means.csv", index=False); RS["T_CUB"] = md(cu_.rename(columns={"conc": "dose"}), "{:.1f}")
zr = pd.read_csv(T / "s11_zinc_rank_consistency.csv"); zr = zr[zr.trait.isin(["a", "lnA"])]; RS["T_ZN_RANK"] = md(zr, "{:.2f}")
un = pd.read_csv(R / "unfiltered/mixed_model_summary.csv"); mm0 = pd.read_csv(R / "mixed_model_summary.csv")
sens = un[["Metal", "model", "dose_effect"]].merge(mm0[["Metal", "model", "dose_effect"]], on=["Metal", "model"], suffixes=(" unfiltered", " main")); sens["model"] = sens.model.map({"M0": "total", "M1": "at fixed size"})
RS["T_SENS_MM"] = md(sens)
U = REP / "sensitivity_unfiltered/tables"; rows = []
om_u = pd.read_csv(U / "s1_omnibus_tests.csv"); om_m = pd.read_csv(T / "s1_omnibus_tests.csv")
for mt in ("Chromium", "Copper", "Iron", "Lead"):
    for mdl in ("baseline_size_adjusted", "population_baseline_size_adjusted"):
        a = om_u[(om_u.Metal == mt) & (om_u.model == mdl)].F.iloc[0]; b = om_m[(om_m.Metal == mt) & (om_m.model == mdl)].F.iloc[0]; rows.append(dict(quantity=f"{mt}: omnibus F, {mdl.replace('_', ' ')}", unfiltered=a, main=b))
s5u = pd.read_csv(U / "s5_split_half.csv"); s5m = pd.read_csv(T / "s5_split_half.csv")
for mt in ("Chromium", "Copper", "Iron", "Lead"):
    a = s5u[(s5u.Metal == mt) & (s5u.measure == "a_adj")].iloc[0]; b = s5m[(s5m.Metal == mt) & (s5m.measure == "a_adj")].iloc[0]
    rows.append(dict(quantity=f"{mt}: split-half gap (size-adjusted)", unfiltered=a.split_gap_mean, main=b.split_gap_mean)); rows.append(dict(quantity=f"{mt}: reliability of strain change", unfiltered=a.reliability_change_A_vs_B, main=b.reliability_change_A_vs_B))
h_u = pd.read_csv(U / "s6_run_heterogeneity.csv"); h_m = pd.read_csv(T / "s6_run_heterogeneity.csv")
for mt in ("Chromium", "Copper", "Iron", "Lead"):
    a = om_u[(om_u.Metal == mt) & (om_u.model == "dose_response_size_adjusted")].F.iloc[0]; b = om_m[(om_m.Metal == mt) & (om_m.model == "dose_response_size_adjusted")].F.iloc[0]; rows.append(dict(quantity=f"{mt}: omnibus F, species by dose", unfiltered=a, main=b))
    a = s5u[(s5u.Metal == mt) & (s5u.measure == "a_adj")].n_strains.iloc[0]; b = s5m[(s5m.Metal == mt) & (s5m.measure == "a_adj")].n_strains.iloc[0]; rows.append(dict(quantity=f"{mt}: strains in the split-half analysis", unfiltered=a, main=b))
p_u = pd.read_csv(U / "s2_phylo_signal.csv"); p_m = pd.read_csv(T / "s2_phylo_signal.csv")
for mt in ("Chromium", "Copper", "Iron", "Lead"):
    f_ = lambda d: d[(d.Metal == mt) & (d.scope == "all_strains_with_tip") & (d.trait == "baseline a*")].pagel_lambda.iloc[0]; rows.append(dict(quantity=f"{mt}: Pagel lambda, baseline a*", unfiltered=f_(p_u), main=f_(p_m)))
RS["T_SENS_OTHER"] = md(pd.DataFrame(rows))
body = open("analysis/carotenoid_stress_vs_baseline/scripts/strat/report_text.md").read()
subs = dict(T_MM=t_mm, T_M2=t_m2, T_COR_ALL=topcor("all_doses"), T_COR_0=topcor("dose0"), T_COR_WITHIN=topcor("within_strain_dose"), T_OM=t_om, T_SB=t_sb, T_SS=t_ss, T_PO=t_po, T_PB=t_pb, T_PP=t_pp,
            T_PH=t_ph, T_SR=t_sr, T_RS=t_rs, T_S4=t_s4, T_S5=t_s5, T_H6=t_h6, T_P6=t_p6, T_V6=t_v6, T_SJ=t_sj)

# ---- section 10 (minimum area)
s8m = pd.read_csv(T / "s8_minarea_models.csv"); s8m = s8m[s8m.model == "M1"]
t_s8m = md(s8m.pivot(index="Metal", columns="min_area", values="dose_effect").reset_index())
s8f = pd.read_csv(T / "s8_minarea_dose_factor.csv"); s8w = pd.read_csv(T / "s8_minarea_wells_by_dose.csv")
def minarea_tab(m):
    e = s8f[s8f.Metal == m].pivot(index="conc", columns="min_area", values="effect_at_fixed_size"); e.columns = [f"effect @{c}" for c in e.columns]
    n = s8w[s8w.Metal == m].pivot(index="conc", columns="min_area", values="wells"); n.columns = [f"wells @{c}" for c in n.columns]
    return md(e.join(n, how="outer").reset_index().rename(columns={"conc": "dose"}), "{:.1f}").replace(".0 |", " |")
s8t = pd.read_csv(T / "s8_minarea_topdose.csv"); t_s8t = md(s8t[["Metal", "min_area", "top_dose", "n_strains_top", "top_dose_mean_a", "top_dose_sd_across_strains", "median_area_top_px"]])
subs.update(T_S8_MODEL=t_s8m, T_S8_CR=minarea_tab("Chromium"), T_S8_PB=minarea_tab("Lead"), T_S8_TOP=t_s8t)
subs.update(RS)
for k, v in subs.items(): body = body.replace("{{" + k + "}}", v)
import re
for f in sorted(set(re.findall(r"\{\{IMG:([a-z0-9_]+)\|([^}]*)\}\}", body))): body = body.replace("{{IMG:%s|%s}}" % f, img(f[0], f[1]))
assert "{{" not in body, re.findall(r"\{\{[^}]*\}\}", body)[:5]
# pandoc needs a blank line before a bullet list that follows a paragraph line (GitHub does not)
_lines = body.split("\n"); _out = []
for _i, _l in enumerate(_lines):
    if _l.startswith("- ") and _out and _out[-1].strip() and not _out[-1].startswith(("- ", "  ", "|")): _out.append("")
    _out.append(_l)
body = "\n".join(_out)
(REP / "REPORT.md").write_text(body); print("REPORT.md written,", len(body.splitlines()), "lines")
