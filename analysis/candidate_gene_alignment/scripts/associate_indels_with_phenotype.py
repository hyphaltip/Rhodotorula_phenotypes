#!/usr/bin/env python3
"""Associate panel-segregating coding INDELs (esp. frameshift/length-changing) at the
candidate genes with color/pigment phenotypes, using the SAME population-aware battery
as associate_variants_with_phenotype.py.

Motivation (2026-08-26 correction): the pilot's per-gene indel screen (screen_indels.sh)
silently reported 0 indels because bcftools wasn't on PATH in that session -- the pipe +
`|| true` swallowed the failure. The real INDEL VCF has hundreds of records in every
candidate gene's CDS+/-2kb, including segrating indels INSIDE the CDS exons of all 10
genes (see results/gene_coding_indel_screen.csv). Those coding indels are absent from
the SNP-only variant tables, so if one of them (e.g. a frameshift in a carotenoid
biosynthesis gene) actually drives a color phenotype, the SNP-only association would
miss it. This script pulls the indel genotypes straight from the INDEL VCF and runs the
same within-population battery.

Usage:
  pixi run python scripts/associate_indels_with_phenotype.py \
    --gene-exons results/gene_coding_indel_screen.json \
    --pop-csv ... --culled-list ... --collapse-work-dir ... --out out.csv

Requires bcftools on PATH (set BCFTOOLS env var to an absolute path to be safe).
"""
import argparse
import gzip
import json
import os
import subprocess
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from associate_variants_with_phenotype import (
    load_pop_map, load_culled_set, load_collapse_maps,
    within_pop_test, meta_analysis, population_covariate_test, veridict,
    query_indel_genotypes,
)

CULLED_DEFAULT = "analysis/gwas/results/gwas/near_clone_culling/culled_keep.txt"
POP_DEFAULT = "analysis/gwas/data/prior_run_state/pop_assignment_at_run.csv"
COLLAPSE_DEFAULT = "analysis/gwas/results/gwas/population_vs_locus/work"
PHENO_DEFAULT = "analysis/gwas/results/gwas_next_phenotypes.csv"
INDEL_VCF = "data/raw/genotypes/RmucY2510_v2/RmucY2510_v2.All.INDEL.combined_selected.vcf.gz"
ACCEPTED = "analysis/gwas/results/strain_reconciliation/accepted_vcf_ids.txt"
GFF = ("/bigdata/stajichlab/shared/projects/Population_Genomics/Rhodotorula_mucilaginosa_NRRLY2510/"
       "genome/Rhodotorula_mucilaginosa_NRRL_Y-2510.gff3.gz")


def load_gene_exons_json(path: str) -> dict[str, list[tuple[str, int, int]]]:
    with open(path) as f:
        data = json.load(f)
    out = {}
    for gene, exons in data.items():
        out[gene] = [(e[0], int(e[1]), int(e[2])) for e in exons]
    return out


def load_all_cds_exons(gff3: str, gene_ids: list[str]) -> dict[str, list[tuple[str, int, int]]]:
    opener = gzip.open if gff3.endswith(".gz") else open
    exons_by_gene: dict[str, list[tuple[str, int, int]]] = {g: [] for g in gene_ids}
    with opener(gff3, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "CDS":
                continue
            for gene in gene_ids:
                if f"Parent={gene}-T1" in fields[8]:
                    exons_by_gene[gene].append((fields[0], int(fields[3]), int(fields[4])))
    for gene in gene_ids:
        assert exons_by_gene[gene], f"no CDS exons for {gene}-T1"
    return exons_by_gene


def iter_segregating_indels(exons_by_gene, accepted, bcftools, per_gene: bool = True):
    """Yield (gene, (pos, ref, alt), {strain: dosage}) for records with >=1 alt in the
    accepted panel."""
    for gene, exons in exons_by_gene.items():
        seen = {}
        for (scf, s, e) in exons:
            records = query_indel_genotypes(INDEL_VCF, scf, s, e, accepted, bcftools=bcftools)
            for key, genos in records.items():
                if key in seen:
                    continue
                seen[key] = True
                n_alt = sum(1 for v in genos.values() if v == "1")
                if n_alt > 0:
                    yield gene, scf, key, genos, n_alt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene-ids", nargs="+", required=True)
    ap.add_argument("--pop-csv", default=POP_DEFAULT)
    ap.add_argument("--culled-list", default=CULLED_DEFAULT)
    ap.add_argument("--collapse-work-dir", default=COLLAPSE_DEFAULT)
    ap.add_argument("--pheno-csv", default=PHENO_DEFAULT)
    ap.add_argument("--traits", nargs="+", default=["chroma", "sat", "bright", "lab_L", "lab_a", "lab_b"])
    ap.add_argument("--bcftools", default=os.environ.get("BCFTOOLS", "bcftools"))
    ap.add_argument("--out", required=True)
    ap.add_argument("--gene-coords", default=None,
                     help="optional JSON {gene: [[scaffold,start,end],...]}; defaults to GFF3 lookup")
    args = ap.parse_args()

    pop_all = load_pop_map(args.pop_csv)
    culled = load_culled_set(args.culled_list)
    rep_maps = load_collapse_maps(args.collapse_work_dir)
    accepted = set(l.strip() for l in open(ACCEPTED) if l.strip())
    pheno = pd.read_csv(args.pheno_csv).set_index("strain_code")

    exons_by_gene = (load_gene_exons_json(args.gene_coords) if args.gene_coords
                     else load_all_cds_exons(GFF, args.gene_ids))

    results = []
    details = []
    for gene, scf, (pos, ref, alt), genos, n_alt_total in iter_segregating_indels(
            exons_by_gene, accepted, args.bcftools):
        ref, alt = str(ref), str(alt)
        # Multiallelic records carry comma-joined ALT; measure per-ALLELE length change
        # (the joined-string length is meaningless). frameshift = ANY allele whose net
        # length change /= 0 mod 3.
        alt_alleles = [a for a in alt.split(",") if a]
        per_allele_len = sorted({len(a) - len(ref) for a in alt_alleles})
        is_frameshift = any((d % 3) != 0 for d in per_allele_len)
        len_desc = "+".join(f"{'+' + str(d) if d >= 0 else d}" for d in per_allele_len)
        kind = ("del" if max(per_allele_len, default=0) < 0 and min(per_allele_len, default=0) <= 0 else
                "ins" if min(per_allele_len, default=0) > 0 else "insdel/snv")
        label = (f"{gene}:{pos} indel {ref}->{alt} "
                 f"({kind} {len_desc}bp, {'frameshift' if is_frameshift else 'in-frame'})")

        for trait in args.traits:
            strain_info = {s: float(1 if gt == "1" else 0) for s, gt in genos.items()
                           if s in culled and s in pop_all and gt in ("0", "1")}
            phen = {s: float(v) for s, v in pheno[trait].items()
                    if pd.notna(v) and s in strain_info}
            testable = set(strain_info) & set(phen)
            if len(testable) < 6:
                results.append(dict(gene=gene, pos=pos, variant=label, ref=ref, alt=alt,
                                    allele_len=",".join(str(d) for d in per_allele_len), frameshift=is_frameshift,
                                    trait=trait, n_strains=len(testable), n_alt=n_alt_total,
                                    meta_beta=np.nan, meta_se=np.nan, meta_p=np.nan,
                                    n_pops_tested=0, n_directionally_consistent=0,
                                    partial_r2=np.nan, partial_p=np.nan, verdict="untested_low_n"))
                continue
            dosage = {s: strain_info[s] for s in testable}
            a_df = within_pop_test(dosage, phen, pop_all, testable,
                                   {p: rep_maps.get(p, {s: s for s in testable})
                                    for p in sorted(set(pop_all.get(s) for s in testable))})
            meta = meta_analysis(a_df)
            cov = population_covariate_test(dosage, phen, pop_all, list(testable))
            verdict = veridict(meta, cov)
            results.append(dict(gene=gene, pos=pos, variant=label, ref=ref, alt=alt,
                                allele_len=",".join(str(d) for d in per_allele_len), frameshift=is_frameshift,
                                trait=trait, n_strains=len(testable), n_alt=n_alt_total,
                                meta_beta=meta["meta_beta"], meta_se=meta["meta_se"], meta_p=meta["meta_p"],
                                n_pops_tested=meta["n_pops"],
                                n_directionally_consistent=meta["n_directionally_consistent"],
                                partial_r2=cov["partial_r2"], partial_p=cov["partial_p"],
                                cov_n=cov["n"], cov_n_pops=cov["n_pops"], verdict=verdict))
            for _, pr in a_df.iterrows():
                details.append(dict(gene=gene, pos=pos, variant=label, trait=trait,
                                    population=pr.population, n_raw=pr.n_raw,
                                    n_effective=pr.n_effective, n0=pr.n0, n1=pr.n1,
                                    beta=pr.beta, se=pr.se, p=pr.p,
                                    tested=pr.tested, reason=pr.reason))

    out = pd.DataFrame(results)
    valid_meta = out["meta_p"].dropna()
    if len(valid_meta) > 0:
        m = len(valid_meta)
        order = valid_meta.sort_values().index
        p_sorted = valid_meta.loc[order].to_numpy()
        ranks = np.arange(1, m + 1)
        fdr = p_sorted * m / ranks
        fdr = np.minimum.accumulate(fdr[::-1])[::-1]
        out["meta_p_fdr"] = pd.Series(fdr, index=order).reindex(out.index).astype(float)
    else:
        out["meta_p_fdr"] = np.nan
    out["fdr_sig"] = out["meta_p_fdr"] < 0.05
    out = out.sort_values(["fdr_sig", "meta_p"], ascending=[False, True]).reset_index(drop=True)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    out.to_csv(args.out, index=False)
    pd.DataFrame(details).to_csv(args.out.replace(".csv", "_within_pop_details.csv"), index=False)

    frameshift = out[out["frameshift"]]
    print(f"\nWrote {args.out} ({len(out)} indel x trait tests)")
    print(f"  frameshift indels tested: {frameshift.shape[0]}; FDR-significant: {out['fdr_sig'].sum()}")
    print("Top 15 by raw meta_p:")
    cols = ["gene", "pos", "variant", "trait", "n_alt", "meta_p", "meta_p_fdr",
            "partial_r2", "n_pops_tested", "n_directionally_consistent", "verdict"]
    print(out.head(15)[cols].to_string(index=False))
    print("\nVerdict counts:", out["verdict"].value_counts().to_dict())


if __name__ == "__main__":
    main()
