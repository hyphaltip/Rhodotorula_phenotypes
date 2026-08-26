#!/usr/bin/env python3
"""Build a per-gene variant table joining snpEff consequence calls with per-strain
genotype, population, phenotype value, and (for GWAS-locus genes) the lead SNP's own
genotype and each coding variant's own Tier A marginal p-value.

Per both reviews (Opus general + fable bioinformatics-expert): this table is
DESCRIPTIVE, not a second association test -- it exists to show, within one gene,
which (if any) coding variant plausibly drives an already-established GWAS signal
versus which are LD-linked passengers. It reports each site's own Tier A p-value for
context, not as a new hypothesis test.

Input: an snpEff-annotated VCF restricted to the gene's CDS+flank region (produced by
running `java -jar $SNPEFFJAR eff ...` on a bcftools-extracted region, see
CANDIDATE_GENE_ALIGNMENT.md), the same VCF's genotypes, and this repo's existing
phenotype/population CSVs.
"""
import argparse
import csv
import gzip
import re
import subprocess

import pandas as pd

ANN_FIELDS = [
    "allele", "consequence", "impact", "gene_name", "gene_id", "feature_type",
    "feature_id", "biotype", "rank", "hgvs_c", "hgvs_p", "cdna_pos", "cds_pos",
    "protein_pos", "distance", "errors",
]


def parse_ann(ann_str: str, target_gene_id: str) -> dict | None:
    """ANN can carry one annotation per overlapping transcript, comma-separated.
    Return the one matching target_gene_id, or None if the site doesn't annotate
    against this gene at all (shouldn't happen given we already restricted to a
    region around it, but guard anyway)."""
    for entry in ann_str.split(","):
        parts = entry.split("|")
        if len(parts) < len(ANN_FIELDS):
            parts = parts + [""] * (len(ANN_FIELDS) - len(parts))
        rec = dict(zip(ANN_FIELDS, parts))
        if rec["gene_id"] == target_gene_id:
            return rec
    return None


def load_assoc_pvalues(assoc_path: str) -> dict[int, float]:
    """GEMMA .assoc.txt: columns include chr, ps (position), p_wald (or p_lrt)."""
    df = pd.read_csv(assoc_path, sep="\t")
    pcol = "p_wald" if "p_wald" in df.columns else ("p_lrt" if "p_lrt" in df.columns else "p_score")
    return dict(zip(df["ps"], df[pcol]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene-id", required=True)
    ap.add_argument("--ann-vcf", required=True, help="snpEff-annotated VCF/VCF.GZ for this gene's region")
    ap.add_argument("--strain-list", required=True)
    ap.add_argument("--pop-csv", required=True, help="Strain,Pop map (e.g. pop_assignment_at_run.csv)")
    ap.add_argument("--pheno-csv", required=True, help="gwas_next_phenotypes.csv (strain_code + trait cols)")
    ap.add_argument("--lead-snp", default=None, help="scaffold:pos of this gene's GWAS lead SNP, or omit for pathway genes with no locus")
    ap.add_argument("--lead-snp-trait", default=None, help="trait name whose assoc file p-values to report")
    ap.add_argument("--assoc-file", default=None, help="Tier A .assoc.txt for --lead-snp-trait, to report each site's own p-value")
    ap.add_argument("--out-csv", required=True)
    args = ap.parse_args()

    strains = sorted(set(l.strip() for l in open(args.strain_list) if l.strip()))
    pop = {}
    with open(args.pop_csv) as f:
        for r in csv.DictReader(f):
            pop[r["Strain"]] = r["Pop"]
    pheno = pd.read_csv(args.pheno_csv).set_index("strain_code")

    pvals = load_assoc_pvalues(args.assoc_file) if args.assoc_file else {}

    lead_scaffold = lead_pos = None
    if args.lead_snp:
        lead_scaffold, lead_pos = args.lead_snp.split(":")
        lead_pos = int(lead_pos)

    opener = gzip.open if args.ann_vcf.endswith(".gz") else open
    rows = []
    lead_snp_genotypes = {}
    with opener(args.ann_vcf, "rt") as f:
        sample_order = None
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                sample_order = line.rstrip("\n").split("\t")[9:]
                continue
            fields = line.rstrip("\n").split("\t")
            chrom, pos, _id, ref, alt, qual, filt, info = fields[:8]
            pos = int(pos)
            gt_fields = fields[9:]
            m = re.search(r"ANN=([^;]+)", info)
            ann = parse_ann(m.group(1), args.gene_id) if m else None

            genos = {}
            for s, g in zip(sample_order, gt_fields):
                gt = g.split(":")[0].replace("|", "/").split("/")[0]
                genos[s] = gt if gt in ("0", "1") else None

            if lead_scaffold and chrom == lead_scaffold and pos == lead_pos:
                lead_snp_genotypes = dict(genos)

            if ann is None:
                continue  # this site doesn't annotate against the target gene (e.g. only hits a neighbor)

            row = {
                "gene_id": args.gene_id,
                "scaffold": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "consequence": ann["consequence"],
                "impact": ann["impact"],
                "hgvs_c": ann["hgvs_c"],
                "hgvs_p": ann["hgvs_p"],
                "tier_a_p_value": pvals.get(pos, None),
                "is_lead_snp": (lead_scaffold == chrom and lead_pos == pos),
            }
            for s in strains:
                gt = genos.get(s)
                row[f"gt__{s}"] = {"0": "ref", "1": "alt", None: "missing"}[gt]
            rows.append(row)

    table = pd.DataFrame(rows)
    table.to_csv(args.out_csv, index=False)
    print(f"{args.gene_id}: {len(table)} annotated CDS-region variant rows -> {args.out_csv}")
    if not table.empty:
        print(f"  consequence counts:\n{table['consequence'].value_counts().to_string()}")

    # Companion per-strain summary: population, phenotype snapshot, lead-SNP genotype
    strain_rows = []
    for s in strains:
        r = {"strain": s, "population": pop.get(s, "NA")}
        if s in pheno.index:
            for trait in ("chroma", "sat", "bright", "lab_L", "lab_a", "lab_b"):
                if trait in pheno.columns:
                    r[trait] = pheno.loc[s, trait]
        if lead_snp_genotypes:
            gt = lead_snp_genotypes.get(s)
            r["lead_snp_genotype"] = {"0": "ref", "1": "alt", None: "missing"}.get(gt, "missing")
        strain_rows.append(r)
    strain_summary_path = args.out_csv.replace(".csv", "_strain_context.csv")
    pd.DataFrame(strain_rows).to_csv(strain_summary_path, index=False)
    print(f"  strain context (population/phenotype/lead-SNP-genotype) -> {strain_summary_path}")


if __name__ == "__main__":
    main()
