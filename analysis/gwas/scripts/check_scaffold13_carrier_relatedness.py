#!/usr/bin/env python3
"""Rule out cryptic within-population relatedness among scaffold_13:810026's 3
carrier strains -- the one gap D-29's carrier-permutation test (GWAS.md section 22)
left open: a tight, closely-related sub-clade sharing both the allele and unrelated
causal variation would look "extreme" in that test too, so a real relatedness check
is needed to distinguish "3 independent lineages sharing an allele" from "a small
clonal/near-clonal group."

Uses the existing genome-wide kinship matrix (GEMMA's centered GRM, already computed
for the Tier A LMM scans) rather than rebuilding a phylogeny -- the GRM already IS a
pairwise relatedness estimate at exactly the resolution needed, and using it keeps
this check consistent with the same relatedness measure the rest of the GWAS
pipeline already trusts.

Usage:
  pixi run python3 analysis/gwas/scripts/check_scaffold13_carrier_relatedness.py \
      --kinship analysis/gwas/results/gwas/grm_conditioning/rebuilt_kinship/kins.cXX.txt \
      --fam analysis/gwas/results/gwas/grm_conditioning/rebuilt_kinship/gwas.pruned.fam \
      --pop-csv analysis/gwas/data/prior_run_state/pop_assignment_at_run.csv \
      --carriers TFCN_17-332C-2 TFCN_17-337P-5 TFCN_86A-12 \
      --out analysis/gwas/results/gwas/rare_variant_validation/scaffold13_carrier_relatedness.csv
"""
from __future__ import annotations

import argparse
import itertools

import numpy as np
import pandas as pd


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--kinship", required=True)
    ap.add_argument("--fam", required=True)
    ap.add_argument("--pop-csv", required=True)
    ap.add_argument("--carriers", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    fam = pd.read_csv(args.fam, sep=r"\s+", header=None,
                       names=["fid", "iid", "pid", "mid", "sex", "pheno"])
    kinship = np.loadtxt(args.kinship)
    # ANALYSIS_OK[runtime-assert]: developer tripwire -- the kinship matrix and
    # .fam file are a matched pair from the same GEMMA kinship build; a shape
    # mismatch means the wrong files were passed, not a data-quality issue.
    assert kinship.shape == (len(fam), len(fam)), (
        f"kinship shape {kinship.shape} does not match {len(fam)} strains in {args.fam}"
    )
    idx = {s: i for i, s in enumerate(fam["iid"])}
    for s in args.carriers:
        # ANALYSIS_OK[runtime-assert]: developer tripwire on a CLI argument
        # (--carriers), checked once at startup, not a data-quality check.
        assert s in idx, f"carrier strain {s!r} not found in {args.fam}"

    popdf = pd.read_csv(args.pop_csv)
    popmap = dict(zip(popdf["Strain"], popdf["Pop"]))
    fam["population"] = fam["iid"].map(popmap)

    rows = []
    for a, b in itertools.combinations(args.carriers, 2):
        pop_a, pop_b = fam.loc[idx[a], "population"], fam.loc[idx[b], "population"]
        same_pop = pop_a == pop_b
        k_ab = float(kinship[idx[a], idx[b]])
        if same_pop:
            pop_idx = fam.index[fam["population"] == pop_a].tolist()
            within = [float(kinship[i, j]) for i, j in itertools.combinations(pop_idx, 2)]
            percentile = float(np.mean(np.array(within) < k_ab) * 100)
            baseline = f"within-population-{pop_a}"
        else:
            offdiag = kinship[np.triu_indices(len(fam), k=1)]
            percentile = float(np.mean(offdiag < k_ab) * 100)
            baseline = "genome-wide"
        rows.append({
            "strain_a": a, "strain_b": b, "pop_a": pop_a, "pop_b": pop_b,
            "same_population": same_pop, "kinship": k_ab,
            "comparison_baseline": baseline, "percentile_in_baseline": percentile,
        })

    out = pd.DataFrame(rows)
    out.to_csv(args.out, index=False)
    print(out.to_string(index=False))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
