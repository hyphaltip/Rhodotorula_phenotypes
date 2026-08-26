#!/usr/bin/env python3
"""Population-confounding sanity check for Tier A top hits, per the quant-genetics
consult (2026-08-25, D-16): for each trait's top FDR-significant SNP, cross-tab
the alt-allele frequency by population (Rmuc_PopAssigned.csv). A locus where one
population carries the alt allele near-fixed while others carry it near-absent
is a population-private-variant pattern -- indistinguishable from a population
main effect rather than a true phenotype association, even after kinship
correction (kinship corrects genome-wide relatedness, not a single locus that
happens to be a population marker).

Flag: population_confound_risk = True if (max per-pop AF - min per-pop AF) > 0.8
AND the SNP is not itself near-fixed genome-wide (0.05 < overall AF < 0.95) --
i.e., a large population AF swing that isn't just a rare/near-fixed variant
everywhere.
"""
import argparse
import csv
import subprocess

import pandas as pd


def genotype_at_site(bcftools_bin: str, vcf: str, scaffold: str, pos: int, samples: list[str]) -> dict[str, float]:
    """Returns {sample: alt_allele_dosage} for a haploid-encoded VCF (0 or 1)."""
    region = f"{scaffold}:{pos}-{pos}"
    out = subprocess.run(
        [bcftools_bin, "query", "-r", region, "-f", "[%SAMPLE=%GT\\t]\\n", vcf],
        check=True, capture_output=True, text=True,
    )
    line = out.stdout.strip()
    assert line, f"no VCF record found at {region}"
    dosage = {}
    for tok in line.split("\t"):
        if not tok or "=" not in tok:
            continue
        s, gt = tok.split("=")
        gt = gt.replace("|", "/").split("/")[0]
        if gt in (".", ""):
            continue
        dosage[s] = float(gt)
    return dosage


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tiera-summary", required=True)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--pop-csv", required=True)
    ap.add_argument("--bcftools", default="bcftools")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    summary = pd.read_csv(args.tiera_summary)
    pop = {}
    with open(args.pop_csv) as f:
        for r in csv.DictReader(f):
            pop[r["Strain"]] = r["Pop"]

    rows = []
    for _, r in summary.iterrows():
        if r["n_fdr05"] == 0:
            continue
        rs = r["top_snp"]
        parts = rs.split(":")
        scaffold, pos = parts[0], int(parts[1])
        dosage = genotype_at_site(args.bcftools, args.vcf, scaffold, pos, [])
        df = pd.DataFrame({"strain": list(dosage.keys()), "dosage": list(dosage.values())})
        df["pop"] = df["strain"].map(pop)
        df = df.dropna(subset=["pop"])
        overall_af = df["dosage"].mean() / 2 if df["dosage"].max() > 1 else df["dosage"].mean()
        per_pop_af = df.groupby("pop")["dosage"].apply(
            lambda s: s.mean() / 2 if s.max() > 1 else s.mean()
        )
        af_range = per_pop_af.max() - per_pop_af.min()
        confound_risk = bool(af_range > 0.8 and 0.05 < overall_af < 0.95)
        rows.append({
            "trait": r["trait"], "top_snp": rs, "n_fdr05": r["n_fdr05"],
            "overall_alt_af": round(overall_af, 3),
            "per_pop_af_min": round(per_pop_af.min(), 3),
            "per_pop_af_max": round(per_pop_af.max(), 3),
            "af_range_across_pops": round(af_range, 3),
            "population_confound_risk": confound_risk,
            "per_pop_af_detail": "; ".join(f"{p}={v:.3f}" for p, v in per_pop_af.items()),
        })
        flag = "FLAGGED" if confound_risk else "ok"
        print(f"{r['trait']}: {rs}  overall_af={overall_af:.3f}  "
              f"per-pop range=[{per_pop_af.min():.3f},{per_pop_af.max():.3f}]  {flag}")

    out = pd.DataFrame(rows)
    out.to_csv(args.out, index=False)
    print(f"\nWrote {args.out}")
    n_flag = out["population_confound_risk"].sum() if len(out) else 0
    print(f"{n_flag}/{len(out)} top hits flagged as population-confound risk")


if __name__ == "__main__":
    main()
