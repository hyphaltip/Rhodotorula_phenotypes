#!/usr/bin/env python3
"""Correlate the new cu_doseauc_v0151 trait against existing copper-response
GWAS traits, for the strain set where both are available. Sanity check run
before -- and reported alongside -- the GEMMA scan for this trait (D-24).

Usage: pixi run python3 analysis/gwas/scripts/compare_copper_v0151_trait.py
"""
from __future__ import annotations

import pathlib

import pandas as pd

REPO = pathlib.Path(__file__).resolve().parents[3]
GWAS = REPO / "analysis/gwas"
PHENO_CSV = GWAS / "results/gwas_next_phenotypes_fam_order.csv"
OUT_CSV = GWAS / "results/gwas/tierA_summary/cu_doseauc_v0151_correlations.csv"

EXISTING_COPPER_TRAITS = [
    "AUC_0", "AUC_10", "AUC_20", "AUC_30", "AUC_ratio_10",
    "resilience_30", "cu_dose_slope", "IC50_est",
]
NEW_TRAIT = "cu_doseauc_v0151"


def main() -> None:
    df = pd.read_csv(PHENO_CSV)
    # ANALYSIS_OK[runtime-assert]: developer tripwire on this pipeline's own
    # fixed-schema intermediate file, checked once at load time.
    assert NEW_TRAIT in df.columns, f"{PHENO_CSV} missing {NEW_TRAIT} -- run add_copper_v0151_trait.py first"
    sub = df[["strain_code", NEW_TRAIT] + EXISTING_COPPER_TRAITS].dropna(subset=[NEW_TRAIT])
    n = len(sub)
    print(f"[compare_v0151] n strains with both {NEW_TRAIT} and fam-order phenotypes = {n}")

    rows = []
    for trait in EXISTING_COPPER_TRAITS:
        # ANALYSIS_OK[sample-filter]: pairwise-complete-case per trait (n_pairs
        # is reported per row) -- IC50_est is undefined for strains that never
        # cross the fitted response threshold, so its n is expected to be much
        # smaller than the other traits'; not a hidden drop.
        pair = sub[[NEW_TRAIT, trait]].dropna()
        rho = pair[NEW_TRAIT].corr(pair[trait], method="spearman")
        r = pair[NEW_TRAIT].corr(pair[trait], method="pearson")
        rows.append({"existing_trait": trait, "n_pairs": len(pair),
                      "spearman_rho": rho, "pearson_r": r})
    out = pd.DataFrame(rows).sort_values("spearman_rho", key=abs, ascending=False)
    out.to_csv(OUT_CSV, index=False)
    print(out.to_string(index=False))
    print(f"\n[compare_v0151] wrote {OUT_CSV}")


if __name__ == "__main__":
    main()
