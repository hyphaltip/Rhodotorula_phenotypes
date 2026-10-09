#!/usr/bin/env python3
"""Orthogonal validation for rare single-SNP GWAS hits that the population-vs-locus
battery (check_population_vs_locus.py, GWAS.md section 13) structurally cannot
validate: alleles too rare to have carriers in >=2 populations can never clear that
battery's required >=2-independent-population replication bar, regardless of whether
the effect is real (D-25, D-27).

Design (GWAS.md section 20 Next Steps): population-stratified exact carrier
permutation. Take the REAL carrier strains' population-membership counts as a fixed
stratum (e.g. "5 carriers, all population 3" or "2 from pop1 + 1 from pop6"). Build
the null by exactly enumerating every possible same-stratified pseudo-carrier set
(the same count drawn from the same population(s)) and computing the same
carriers-mean-phenotype statistic for each. Report the observed carriers' percentile
within that exact null distribution.

This differs from the existing battery in a load-bearing way: it does NOT require
independent replication in >=2 populations (impossible for a variant this rare) --
it asks whether the observed carriers are unusual GIVEN exactly where they sit in
the population structure, which is answerable even from a single population. It does
NOT resolve finer-scale relatedness within a population (a tight sub-clade sharing
both the allele and unrelated causal variation would still look "extreme" here) --
that caveat is inherent to any method with this few carriers and is reported
alongside the result, not hidden.

Carrier counts here are always small enough (<=5 total) to enumerate the exact null
exhaustively -- no Monte Carlo, no simulation noise.

Usage:
  pixi run python3 analysis/gwas/scripts/rare_variant_carrier_permutation_test.py \
      --vcf <path> --trait cu_doseauc_v0151 --snp scaffold_9:704260 \
      --pheno-csv analysis/gwas/results/gwas_next_phenotypes_fam_order.csv \
      --pop-csv analysis/gwas/data/prior_run_state/pop_assignment_at_run.csv \
      --fam analysis/gwas/results/gwas/grm_conditioning/rebuilt_kinship/gwas.fam \
      --bcftools <path> --out <path>
"""
from __future__ import annotations

import argparse
import itertools
import subprocess

import numpy as np
import pandas as pd

MAX_ENUMERATE = 5_000_000  # hard safety cap; both loci tested here are well under this


def carriers_for_snp(bcftools: str, vcf: str, snp: str, fam_ids: set[str]) -> pd.DataFrame:
    scaffold, pos = snp.split(":")
    region = f"{scaffold}:{pos}-{pos}"
    out = subprocess.run([bcftools, "query", "-r", region, "-f", "[%SAMPLE=%GT\n]", vcf],
                          check=True, capture_output=True, text=True)
    rows = []
    for line in out.stdout.splitlines():
        if "=" not in line:
            continue
        s, gt = line.split("=", 1)
        if s in fam_ids:
            rows.append({"strain": s, "gt": gt})
    df = pd.DataFrame(rows)
    # ANALYSIS_OK[runtime-assert]: developer tripwire -- the fam panel is a fixed
    # 213-strain list and every member must have a genotype call at this site.
    assert len(df) == len(fam_ids), f"{snp}: expected {len(fam_ids)} genotyped strains, got {len(df)}"
    df["is_carrier"] = df["gt"].isin({"1", "1/1", "1|1"})
    return df


def exact_null(strat_pool: dict[str, np.ndarray], strat_k: dict[str, int]) -> np.ndarray:
    """Enumerate every possible pseudo-carrier draw matching the observed
    per-population carrier counts exactly, return the array of mean-phenotype
    statistics for all such draws."""
    per_pop_sums = []
    per_pop_counts = []
    for pop, k in strat_k.items():
        vals = strat_pool[pop]
        n_combos = len(vals) if k == 1 else int(np.prod([len(vals) - i for i in range(k)]) //
                                                  np.prod(range(1, k + 1)))
        assert n_combos <= MAX_ENUMERATE, (
            f"pop {pop}: choose({len(vals)},{k}) = {n_combos} exceeds enumeration cap"
        )
        sums = np.array([sum(c) for c in itertools.combinations(vals, k)], dtype=float)
        per_pop_sums.append(sums)
        per_pop_counts.append(k)

    total = per_pop_sums[0]
    for s in per_pop_sums[1:]:
        # outer-sum across independent population strata -> cartesian product of draws
        total = (total[:, None] + s[None, :]).ravel()
    n_total_carriers = sum(per_pop_counts)
    return total / n_total_carriers


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--trait", required=True)
    ap.add_argument("--snp", required=True, help="scaffold:pos, e.g. scaffold_9:704260")
    ap.add_argument("--pheno-csv", required=True)
    ap.add_argument("--pop-csv", required=True)
    ap.add_argument("--fam", required=True)
    ap.add_argument("--bcftools", default="bcftools")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    fam = pd.read_csv(args.fam, sep=r"\s+", header=None,
                       names=["fid", "iid", "pid", "mid", "sex", "pheno"])
    popdf = pd.read_csv(args.pop_csv)
    popmap = dict(zip(popdf["Strain"], popdf["Pop"]))
    pheno_df = pd.read_csv(args.pheno_csv)
    pheno_map = dict(zip(pheno_df["strain_code"], pheno_df[args.trait]))

    g = carriers_for_snp(args.bcftools, args.vcf, args.snp, set(fam["iid"]))
    g["population"] = g["strain"].map(popmap)
    g["pheno"] = g["strain"].map(pheno_map)
    n_missing_pheno = g["pheno"].isna().sum()
    print(f"[{args.trait} @ {args.snp}] {len(g)} strains, {g['is_carrier'].sum()} carriers, "
          f"{n_missing_pheno} missing phenotype (excluded)")
    # ANALYSIS_OK[sample-filter]: drops strains with no phenotype value or no
    # population assignment -- count logged above (n_missing_pheno); these
    # strains cannot contribute to either the observed statistic or the
    # exact-enumeration null pool, so they must be excluded, not imputed.
    g = g.dropna(subset=["pheno", "population"])

    carriers = g[g["is_carrier"]]
    strat_k = carriers["population"].value_counts().to_dict()
    print(f"  carrier population strata: {strat_k}")

    strat_pool = {pop: g.loc[g["population"] == pop, "pheno"].to_numpy() for pop in strat_k}
    for pop, k in strat_k.items():
        assert len(strat_pool[pop]) >= k, f"population {pop} has fewer strains than carriers to draw"

    observed_mean = carriers["pheno"].mean()
    null_means = exact_null(strat_pool, strat_k)
    n_null = len(null_means)

    # two-sided: how extreme is the observed mean vs. the exact null distribution,
    # centered on the null's own mean (not on 0 -- the statistic is a raw phenotype
    # mean, whose null center is population-baseline-dependent by construction).
    null_center = null_means.mean()
    obs_dev = abs(observed_mean - null_center)
    null_dev = np.abs(null_means - null_center)
    p_two_sided = (np.sum(null_dev >= obs_dev) + 1) / (n_null + 1)
    percentile = float((null_means < observed_mean).mean() * 100)

    result = pd.DataFrame([{
        "trait": args.trait, "snp": args.snp,
        "n_carriers": int(carriers.shape[0]), "carrier_strata": str(strat_k),
        "carrier_strains": ";".join(carriers["strain"]),
        "observed_carrier_mean": observed_mean,
        "null_mean": null_center, "null_sd": float(null_means.std()),
        "n_null_enumerated": n_null,
        "observed_percentile_in_null": percentile,
        "p_two_sided_exact": p_two_sided,
    }])
    result.to_csv(args.out, index=False)
    print(result.to_string(index=False))
    print(f"  wrote {args.out}")


if __name__ == "__main__":
    main()
