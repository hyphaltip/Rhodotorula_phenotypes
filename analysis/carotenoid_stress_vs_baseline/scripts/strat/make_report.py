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
for f in ("mixed_model_summary", "dose_factor_effects_M2", "wells_lost_by_dose", "astar_trait_correlations", "species_tests"):
    shutil.copy(R / (f + ".csv"), T / (f + ".csv"))
def md(df, floatfmt="{:.2f}"):
    df = df.copy(); df.columns = [str(c) for c in df.columns]; cols = list(df.columns)
    for c in cols:
        if df[c].dtype.kind == "f": df[c] = df[c].map(lambda v: "" if pd.isna(v) else ("<0.001" if (c.lower().startswith("p") and 0 <= v < 0.001) else floatfmt.format(v)))
    rows = [[str(x) for x in r] for r in df.itertuples(index=False)]
    wd = [min(40, max(len(c), *(len(r[i]) for r in rows))) if rows else len(c) for i, c in enumerate(cols)]  # proportional separators -> pandoc column widths
    out = ["| " + " | ".join(cols) + " |", "|" + "|".join("-" * max(3, w) for w in wd) + "|"]
    out += ["| " + " | ".join(r) + " |" for r in rows]
    return "\n".join(out)
img = lambda f, cap: "![](figures/%s.png)\n\n*%s*\n" % (f, cap.replace("\\*", "*").replace("*", "\\*"))   # empty alt text: pandoc would otherwise add a second caption in the PDF
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
    p = cor[cor.scope == sc].pivot(index="trait", columns="Metal", values="rho"); p["mean |rho|"] = p.abs().mean(axis=1); return md(p.sort_values("mean |rho|", ascending=False).head(n).reset_index())
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
s4 = pd.read_csv(T / "s4_size_matched.csv"); t_s4 = md(s4[["Metal", "size_bin", "n_wells", "n_doses", "mean_a", "sd_a", "slope_full_range", "se", "p"]])
s5 = pd.read_csv(T / "s5_split_half.csv"); s5 = s5[s5.measure == "a_adj"]
t_s5 = md(s5[["Metal", "top_dose", "n_strains", "naive_gap", "naive_ci_lo", "naive_ci_hi", "split_gap_mean", "split_gap_ci_lo", "split_gap_ci_hi", "reliability_change_A_vs_B", "reliability_baseline_A_vs_B"]])
h6 = pd.read_csv(T / "s6_run_heterogeneity.csv"); t_h6 = md(h6); p6 = pd.read_csv(T / "s6_per_run_dose_effect.csv"); t_p6 = md(p6); v6 = pd.read_csv(T / "s6_variance_components.csv"); t_v6 = md(v6[["Metal", "n_runs", "share_strain", "share_run", "share_plate", "share_resid", "singular"]])
sj = pd.read_csv(T / "s2_jitter_sensitivity.csv"); t_sj = md(sj[["Metal", "jitter", "cond", "lam", "K", "pK"]], "{:.3g}")
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
mt = open("analysis/carotenoid_stress_vs_baseline/scripts/strat/minarea_text.md").read()
body = body.replace("{{MINAREA}}", mt)
subs.update(T_S8_MODEL=t_s8m, T_S8_CR=minarea_tab("Chromium"), T_S8_PB=minarea_tab("Lead"), T_S8_TOP=t_s8t)
for k, v in subs.items(): body = body.replace("{{" + k + "}}", v)
import re
for f in sorted(set(re.findall(r"\{\{IMG:([a-z0-9_]+)\|([^}]*)\}\}", body))): body = body.replace("{{IMG:%s|%s}}" % f, img(f[0], f[1]))
assert "{{" not in body, re.findall(r"\{\{[^}]*\}\}", body)[:5]
(REP / "REPORT.md").write_text(body); print("REPORT.md written,", len(body.splitlines()), "lines")
