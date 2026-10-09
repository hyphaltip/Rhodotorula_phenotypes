#!/usr/bin/env python3
"""Reconcile and merge the copper-heavy-metal-screen-v0.15.1 dose-response AUC
trait into the existing GWAS phenotype table (fam order).

Strain identity join, not fuzzy: the new dataset's `Strain` column (e.g.
TFCN_17-291Y-1) and this project's `strain_code` are the *same* ID space --
both trace back to Copper.Strain_info.csv's `Strain` column (data/metadata/
Copper.Strain_info.csv), which is also the source of the existing
copper-colony-measurements dataset's SAMPLE_NAME field. strain_code itself is
the post-VCF-reconciliation accepted ID (analysis/gwas/results/
strain_reconciliation/strain_match_table.reviewed.csv), which for the mucilaginosa
panel is an identity map (phenotype_strain_id == vcf_sample_id, TFCN_*) -- see
D-24. No fuzzy matching needed; this script performs and REPORTS an exact-string
join and asserts there are no duplicate claims, so a real mismatch (e.g. an ID
typo) fails loudly rather than silently dropping strains.

Output: adds one new column (cu_doseauc_v0151) to
analysis/gwas/results/gwas_next_phenotypes_fam_order.csv (213-strain fam order,
already used by run_tiera_gemma.sh), and writes a join report to
analysis/gwas/results/strain_reconciliation/copper_v0151_match_table.csv.

Usage: pixi run python3 analysis/gwas/scripts/add_copper_v0151_trait.py
"""
from __future__ import annotations

import pathlib

import pandas as pd

REPO = pathlib.Path(__file__).resolve().parents[3]
GWAS = REPO / "analysis/gwas"
FAM_PHENO = GWAS / "results/gwas_next_phenotypes_fam_order.csv"
NEW_TRAIT_SRC = REPO / "data/raw/copper-heavy-metal-screen-v0.15.1/copper_auc_mean_by_strain.csv"
MATCH_OUT = GWAS / "results/strain_reconciliation/copper_v0151_match_table.csv"
TRAIT_COL = "cu_doseauc_v0151"


def main() -> None:
    fam_pheno = pd.read_csv(FAM_PHENO)
    # ANALYSIS_OK[runtime-assert]: developer tripwire on a fixed pipeline-internal
    # file's schema, not user input; this script is never run under python -O.
    assert "strain_code" in fam_pheno.columns, f"{FAM_PHENO} missing strain_code column"
    n_fam = len(fam_pheno)
    print(f"[add_copper_v0151] fam-order phenotype table: {n_fam} strains")

    new = pd.read_csv(NEW_TRAIT_SRC)
    # ANALYSIS_OK[runtime-assert]: same rationale as above -- fixed ingested
    # dataset's schema/join-key uniqueness, checked once at load time.
    for col in ("Strain ID", "Strain", "mean_auc_rate"):
        assert col in new.columns, f"{NEW_TRAIT_SRC} missing expected column {col!r}"
    n_new = len(new)
    n_dup_strain = new["Strain"].duplicated().sum()
    # ANALYSIS_OK[runtime-assert]: developer tripwire on the ingested dataset's
    # join-key uniqueness, not user input; script never run under python -O.
    assert n_dup_strain == 0, (
        f"{NEW_TRAIT_SRC} has {n_dup_strain} duplicate 'Strain' values -- "
        "join key must be unique before merging"
    )
    print(f"[add_copper_v0151] source trait table: {n_new} strains, "
          f"mean_auc_rate range [{new.mean_auc_rate.min():.3f}, {new.mean_auc_rate.max():.3f}]")

    merged = fam_pheno[["strain_code"]].merge(
        new[["Strain ID", "Strain", "mean_auc_rate"]],
        left_on="strain_code", right_on="Strain", how="left",
        validate="one_to_one",
    )
    n_matched = merged["mean_auc_rate"].notna().sum()
    n_unmatched = n_fam - n_matched
    print(f"[add_copper_v0151] join: {n_matched}/{n_fam} fam-order strains matched "
          f"({n_unmatched} unmatched -- no copper_v0.15.1 AUC for these strains)")

    match_report = merged.rename(columns={"mean_auc_rate": "cu_doseauc_v0151"}).copy()
    match_report["matched"] = match_report["cu_doseauc_v0151"].notna()
    match_report.to_csv(MATCH_OUT, index=False)
    print(f"[add_copper_v0151] wrote join report: {MATCH_OUT}")

    fam_pheno[TRAIT_COL] = merged["mean_auc_rate"].to_numpy()
    fam_pheno.to_csv(FAM_PHENO, index=False)
    print(f"[add_copper_v0151] wrote updated phenotype table "
          f"(added column {TRAIT_COL}): {FAM_PHENO}")


if __name__ == "__main__":
    main()
