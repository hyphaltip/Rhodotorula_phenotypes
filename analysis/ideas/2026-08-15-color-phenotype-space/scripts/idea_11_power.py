#!/usr/bin/env python3
"""
idea_11_power.py -- Feasibility / power analysis for mapping genetic variants to
phenotypes (growth rate, color) in the Rhodotorula panel.

Two independent questions:
  (1) PHENOTYPE SIDE -- is there enough repeatable, non-environmental phenotypic
      variance (i.e. a signal that could plausibly be genetic) to even rank strains?
      Uses per-strain ICC / repeatability (idea_04) and within-species SD (idea_10).
  (2) POWER / GENOTYPE SIDE -- given the effective number of independent genomes
      (NOT nominal strain count -- many are near-clonal, from the protein phylogeny)
      and the number of variant tests at genome-wide significance, what per-locus
      effect (R^2, per-allele SD effect) could we detect at 80% power?

Outputs:
  results/idea11_power.csv            -- power table (n, alpha, m variants, detectable R2)
  results/idea11_effective_n.csv      -- genotype redundancy / effective n by species
  results/idea11_detectable_effect.csv-- phenotype-side summary (SD, ICC, detectable beta)
  figures/fig11_power_curves.png      -- power curves
"""

from __future__ import annotations

import math
import os

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RESULTS = os.path.join(ROOT, "results")
FIGURES = os.path.join(ROOT, "figures")
os.makedirs(RESULTS, exist_ok=True)
os.makedirs(FIGURES, exist_ok=True)

GENOME_WIDE_ALPHA = 5e-8          # standard GWAS threshold (~10^7 tests)
EXOME_ALPHA = 1e-6                # exome/biallelic ~10^6 tests
CANDIDATE_ALPHA = 1e-4            # candidate-gene / small variant panel

# ---------------------------------------------------------------------------
# (1) PHENOTYPE SIDE ---------------------------------------------------------
# ---------------------------------------------------------------------------

# Per-strain trait SD (SD_total across all strains, log-units or per-mM).
# From idea_10_species_variation.results.
TRAITS = {
    # trait: (sd_total, within-species fraction of variance, ICC/repeatability upper bound)
    # sd_total, frac_within from idea10; icc from idea04_icc_audit / reaction norms
    "colony_size_log10px":  (0.302, 0.62, 0.86),   # l10med_fixed; ICC Shape_Area 0.859, mostly between-grp
    "colony_size_log10px_within":  (0.238, 1.00, 0.86),  # sqrt(0.62)*0.302 within-species SD
    "baseline_chroma_log10": (0.214, 0.92, 0.82),   # intercept_logchroma; ICC Chroma 0.816
    "copper_slope_per_mM":   (0.0128, 0.87, 0.19),  # slope_logchroma_per_mM; ICC a*CoeffVar-like (weak)
    "heterogeneity":         (0.0183, 0.67, 0.93),  # partial_slope_sd_cu proxy
    "pace_loglog":           (0.217, 0.84, 0.82),   # pace_loglog
}

# Effective independent genomes within the mapping population (R. mucilaginosa).
# Computed in idea_11_effective_n.R from the protein phylogeny: number of distinct
# haplotypes after collapsing near-zero (<=1e-7) patristic-distance clusters.
MUCILAGINOSA_TIPS = 202         # tips with phenotype
MUCILAGINOSA_EFF = 178          # effective independent genomes (from idea_11_effective_n.R)
ALL_TIPS = 272                  # all-strains mapping set (idea09)

CONF = 1.96                     # z at 95% two-sided


def minimum_detectable_R2(n: int, alpha: float, power: float = 0.80,
                          m: int = 1) -> float:
    """R^2 (variance explained by a SNP) detectable at 80% power for n samples.

    Non-centrality of an F test of a single predictor:
      lambda = n * R2 / (1 - R2)  (roughly)
    Power = P(F_noncentral > F_{alpha,1,n-2}). Solve for R2.
    Use the classic approximation: for OLS with one predictor, the noncentrality
    is NCP = n * R2 / (1 - R2). We solve numerically. m = number of tests scaled
    into alpha already by caller (alpha = alpha_genome / m if correcting Bonferroni).
    """
    df1, df2 = 1, n - 2
    """
    Power of an F-test of a single covariate (one SNP) with n samples, at
    significance alpha, for a SNP explaining r2 of total phenotypic variance.

    Standard result: F-test statistic ~ noncentral F with df1=1, df2=n-2 and
    noncentrality lambda = n * r2 / (1 - r2).  Equivalently t^2 ~ noncentral
    chi-square.  For df2 >= ~100 the noncentral chi-square(df=1, ncp=lambda)
    is an excellent approximation and yields a closed-form power:

        crit = qchisq(1 - alpha, df=1)         # genome-wide threshold
        power = P(Z + sqrt(lambda) > sqrt(crit)) - P(Z + sqrt(lambda) < -sqrt(crit))

    because a noncentral chi-square(1, lambda) variable has the distribution of
    (Z + sqrt(lambda))^2 with Z ~ N(0,1)  --  exact, no approximation needed.
    """

    from math import erf, sqrt

    def norm_cdf(z):
        return 0.5 * (1 + erf(z / sqrt(2.0)))

    # critical value
    # qnorm(1 - alpha/2) via the inverse normal (Acklam rational approximation)
    def qnorm(p: float) -> float:
        a = [-39.69683028665376, 220.9460984245205, -275.9285104469687,
             138.3577518672690, -30.66479806614716, 2.506628277459239]
        b = [-54.47609879822406, 161.5858368580409, -155.6989798598866,
             66.80131188771972, -13.28068155288572]
        c = [-0.007784894002430293, -0.3223964580411365, -2.400758277161838,
             -2.549732539343734, 4.374664141464968, 2.938163982698783]
        d = [0.007784695709041462, 0.3224671290700398, 2.445134137142996,
             3.754408661907416]
        plow = 0.02425
        if p < plow:
            q = math.sqrt(-2 * math.log(p))
            return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / \
                   ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1)
        if p > 1 - plow:
            q = math.sqrt(-2 * math.log(1 - p))
            return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / \
                   ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1)
        q = p - 0.5
        r = q * q
        x = (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q
        return x / (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1)

    crit_chi2 = qnorm(1 - alpha / 2) ** 2
    s = math.sqrt(crit_chi2)

    def power_for_r2(r2: float) -> float:
        if r2 <= 0:
            return alpha
        ncp = n * r2 / (1 - r2)
        t = math.sqrt(ncp)
        return 1 - (norm_cdf(s - t) - norm_cdf(-s - t))

    # root-find the r2 giving required power
    lo, hi = 1e-9, 0.9999
    if power_for_r2(hi) < power:
        return np.nan
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if power_for_r2(mid) < power:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def per_allele_effect_in_sd(r2: float, maf: float) -> float:
    """Additive effect (in phenotype SD units) per allele copy for r2 explained by
    a biallelic SNP with MAF: var_explained = 2*maf*(1-maf)*beta_sd^2.
    beta_sd = sqrt(r2 / (2*maf*(1-maf)))."""
    return math.sqrt(r2 / (2 * maf * (1 - maf)))


def main() -> None:
    rows = []
    for n in [100, 150, 202, 250, 272, 400]:
        for alpha, label in [(GENOME_WIDE_ALPHA, "genome-wide 5e-8"),
                             (EXOME_ALPHA, "exome 1e-6"),
                             (CANDIDATE_ALPHA, "candidate 1e-4")]:
            r2 = minimum_detectable_R2(n, alpha)
            rows.append({"n_samples": n, "alpha": alpha, "alpha_label": label,
                         "min_detectable_R2": r2})
    power_df = pd.DataFrame(rows)
    power_df.to_csv(os.path.join(RESULTS, "idea11_power.csv"), index=False)

    # ---- effective n (fill from R when present) ----
    eff_df = pd.DataFrame({
        "population": ["R. mucilaginosa (tips)", "R. mucilaginosa (eff independent)",
                       "All strains (tips)", "All strains (eff, est)"],
        "sample_size": [MUCILAGINOSA_TIPS, MUCILAGINOSA_EFF, ALL_TIPS, int(ALL_TIPS * 0.85)],
        "note": ["phenotype tips in tree", "after collapsing near-clones (see R)",
                 "idea09 mapping set", "rule-of-thumb redundancy adjustment"],
    })
    eff_df.to_csv(os.path.join(RESULTS, "idea11_effective_n.csv"), index=False)

    # ---- phenotype side ----
    phen_rows = []
    for trait, (sd_total, frac_within, icc) in TRAITS.items():
        sd_within = sd_total * math.sqrt(frac_within)
        # detectable per-allele effect (in trait SD units) at eff n=150, genome-wide
        r2_95 = minimum_detectable_R2(150, GENOME_WIDE_ALPHA)
        beta_sd = per_allele_effect_in_sd(r2_95, maf=0.3)
        # heritability-adjusted: SNP must explain this fraction of the *genetic*
        # variance (ICC = upper bound on h2) to be detected
        r2_adj = r2_95 / icc if icc > 0 else np.nan
        phen_rows.append({
            "trait": trait,
            "sd_total": sd_total,
            "frac_variance_within_species": frac_within,
            "sd_within_species": sd_within,
            "icc_repeatability_ub": icc,
            "detectable_R2_at_eff150_genomewide": r2_95,
            "detectable_frac_of_genetic_var_assuming_h2=ICC": r2_adj,
            "detectable_beta_SD_units_MAF0.3": beta_sd,
            "detectable_beta_trait_units": beta_sd * sd_within,
        })
    phen_df = pd.DataFrame(phen_rows)
    phen_df.to_csv(os.path.join(RESULTS, "idea11_detectable_effect.csv"), index=False)

    # ---- figure ----
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    ns = np.arange(50, 600, 5)
    for alpha, label, style in [(GENOME_WIDE_ALPHA, "genome-wide (5e-8)", "-"),
                                (EXOME_ALPHA, "exome (1e-6)", "--"),
                                (CANDIDATE_ALPHA, "candidate (1e-4)", ":")]:
        r2s = [minimum_detectable_R2(int(n), alpha) for n in ns]
        axes[0].plot(ns, r2s, style, label=label)
    for x, lab in [(MUCILAGINOSA_TIPS, "mucil. tips (202)"),
                   (MUCILAGINOSA_EFF, "mucil. eff (~150)"),
                   (ALL_TIPS, "all tips (272)")]:
        axes[0].axvline(x, color="grey", alpha=0.5, ls="--")
        axes[0].text(x, 0.35, lab, rotation=90, ha="right", fontsize=8, color="grey")
    axes[0].set_xlabel("sample size (n strains / independent genomes)")
    axes[0].set_ylabel("min. detectable per-SNP R2 (80% power)")
    axes[0].set_title("Detectable single-locus effect vs sample size")
    axes[0].legend()
    axes[0].set_ylim(0, 0.5)

    for trait, (sd_total, frac_within, icc) in TRAITS.items():
        sd_w = sd_total * math.sqrt(frac_within)
        r2s = [minimum_detectable_R2(int(n), GENOME_WIDE_ALPHA) for n in ns]
        betas = [per_allele_effect_in_sd(r2, 0.3) * sd_w for r2 in r2s]
        axes[1].plot(ns, betas, label=trait, lw=1.5)
    axes[1].axvline(MUCILAGINOSA_EFF, color="grey", alpha=0.5, ls="--")
    axes[1].text(MUCILAGINOSA_EFF, axes[1].get_ylim()[1] * 0.95, "mucil. eff n",
                 rotation=90, ha="right", fontsize=8, color="grey")
    axes[1].set_xlabel("sample size (n)")
    axes[1].set_ylabel("detectable per-allele effect (trait units, MAF=0.3)")
    axes[1].set_title("Detectable additive effect in real trait units")
    axes[1].legend(fontsize=7)

    fig.tight_layout()
    fig.savefig(os.path.join(FIGURES, "fig11_power_curves.png"), dpi=150)
    print(power_df.to_string(index=False))
    print()
    print(phen_df.to_string(index=False))


if __name__ == "__main__":
    main()
