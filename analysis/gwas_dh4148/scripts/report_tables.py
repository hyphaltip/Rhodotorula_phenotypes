# executed by 12_make_report.py; defines TAB (placeholder -> markdown table)
TAB = {}
r = lambda f: pd.read_csv(T / f)
vc = r("dose0_variance_components.csv"); vc = vc[vc.scope.isin(["Cr+Cu+Pb pooled", "Chromium", "Copper", "Lead", "Cr+Cu+Pb pooled, GWAS panel strains"])]
TAB["T_D0_VC"] = md(vc[["scope", "trait", "n_wells", "n_strains", "share_strain", "share_metal", "share_run", "share_plate", "share_resid"]])
TAB["T_D0_COR"] = md(r("dose0_between_screen_correlations.csv")[["trait", "screen_1", "screen_2", "n_strains", "spearman", "pearson"]])
h = r("dose0_strain_heterogeneity.csv").head(12); h["species"] = h.species.str.replace("Rhodotorula ", "R. "); TAB["T_D0_HET"] = md(h[["sample_name", "species", "n_wells", "a_mean", "a_sd", "a_sd_z"]])
TAB["T_RGR"] = md(r("relative_growth_rate_by_dose.csv"))
s = r("ic50_status_by_metal.csv"); s.columns = [c.split(" (")[0] for c in s.columns]; TAB["T_IC50_STATUS"] = md(s)
g = r("ic50_by_group.csv"); TAB["T_IC50_GRP"] = md(g)
TAB["T_PHENO"] = md(r("gwas_trait_summary.csv"), "{:.3f}")
TAB["T_LIN_THR"] = md(r("lineages_by_threshold.csv")); TAB["T_LIN_SIZES"] = md(r("lineage_sizes.csv"))
ic = r("lineage_icc.csv"); piv = ic.pivot(index="trait", columns="version", values="icc_lineage").reset_index(); pa = ic[ic.version == "run-adjusted"].set_index("trait").anova_p
piv["p (run-adjusted, ANOVA)"] = piv.trait.map(pa).values; piv.columns = ["trait", "ICC raw", "ICC run-adjusted", "p (run-adjusted, ANOVA)"]; TAB["T_ICC"] = md(piv, "{:.2f}")
inf = {}
for tag, lab in (("_unadjusted", "LMM (K only)"), ("_runadj", "LMM, run-adjusted"), ("_runadj_lineage", "LMM, run-adjusted + lineage")):
    d = r(f"gwas_inflation_pve{tag}.csv").set_index("trait"); inf[lab] = d.lambda_lmm; inf["OLS, no relatedness"] = d.lambda_naive_ols
    inf[f"Bonferroni SNPs: {lab}"] = d.n_bonferroni
I = pd.DataFrame(inf).reset_index(); TAB["T_INFL"] = md(I[["trait", "OLS, no relatedness", "LMM (K only)", "LMM, run-adjusted", "LMM, run-adjusted + lineage", "Bonferroni SNPs: LMM (K only)", "Bonferroni SNPs: LMM, run-adjusted", "Bonferroni SNPs: LMM, run-adjusted + lineage"]], "{:.2f}").replace(".00 |", " |")
TAB["T_LINCHECK"] = md(r("gwas_lead_snp_lineage_check.csv"))
l = r("gwas_loci_runadj_lineage.csv"); l["product"] = l["product"].fillna("").str.slice(0, 50)
TAB["T_LOCI_LIN"] = md(l[["trait", "chr", "pos", "af", "beta", "p_wald", "effect", "gene_id", "product", "distance_bp"]].assign(pos=lambda x: x.pos.astype(int)), "{:.3g}")
v = r("gwas_variance_components_additive_epistatic.csv"); TAB["T_VC_EPI"] = md(v[["trait", "n", "h2_additive_only", "h2_additive", "h2_epistatic_AA", "residual_share", "LRT", "p_epistatic"]], "{:.3f}")
e = r("gwas_epistasis_summary.csv"); TAB["T_EPI"] = md(e[["trait", "n_pairs_tested", "bonferroni_threshold", "min_p", "lambda_gc", "n_bonferroni", "n_expected_at_1e_5", "n_at_1e_5"]], "{:.3g}")
u = r("gwas_loci_unadjusted.csv"); b = u[u.bonferroni_sig].copy(); b["product"] = b["product"].fillna("").str.slice(0, 40)
TAB["T_LOCI_UNADJ_TOP"] = md(b.sort_values("p_wald").head(12)[["trait", "chr", "pos", "af", "beta", "p_wald", "n_suggestive_snps_in_locus", "lineages_carrying_alt", "lineage_level_marker", "gene_id", "product"]].assign(pos=lambda x: x.pos.astype(int)), "{:.3g}")
