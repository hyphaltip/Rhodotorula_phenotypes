#!/usr/bin/env python3
"""Screen each candidate gene's actual CDS exons against the INDEL VCF, restricted to
the 213-strain accepted panel, and report whether ANY sung segregating coding indel
invalidates the SNP-only substitution premise for that gene.

The per-gene `indel_screen.csv` from screen_indels.sh counts ALL indels in CDS +/- 2kb
flank -- including intergenic and overlapping another gene's exons. That overstates the
risk. This script computes the precise count: indels whose POS falls inside one of THIS
gene's CDS exons AND that have >=1 alternate call among the 213 accepted strains. A
panel-segregating coding indel means the SNP-only CDS/protein/variant-table for that
gene is an approximation (the sequence is built length-fixed from the SNP VCF; the
indel carriers genuinely have a different-length CDS that this pipeline does not build --
exactly the blind spot the Opus review flagged, which the pilot resolved only by
checking it held 0).

Usage: pixi run python scripts/check_gene_coding_indels.py
       (needs bcftools+hmm on PATH or BCFTOOLS env var)
"""
import argparse
import gzip
import os
import subprocess
import sys

BCFTOOLS = os.environ.get("BCFTOOLS", "bcftools")

GFF = ("/bigdata/stajichlab/shared/projects/Population_Genomics/Rhodotorula_mucilaginosa_NRRLY2510/"
       "genome/Rhodotorula_mucilaginosa_NRRL_Y-2510.gff3.gz")
INDEL_VCF = "data/raw/genotypes/RmucY2510_v2/RmucY2510_v2.All.INDEL.combined_selected.vcf.gz"
ACCEPTED = "analysis/gwas/results/strain_reconciliation/accepted_vcf_ids.txt"

GENES = {
    "OM429_003333": (),
    "OM429_003336": (),
    "OM429_001533": (),
    "OM429_001415": (),
    "OM429_001430": (),
    "OM429_003729": (),
    "OM429_002663": (),
    "OM429_005034": (),
    "OM429_001521": (),
    "OM429_000065": (),
}


def load_cds_exons(gff3: str, gene_id: str) -> list[tuple[str, int, int]]:
    opener = gzip.open if gff3.endswith(".gz") else open
    exons = []
    with opener(gff3, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "CDS":
                continue
            if f"Parent={gene_id}-T1" not in fields[8]:
                continue
            exons.append((fields[0], int(fields[3]), int(fields[4])))
    assert exons, f"no CDS rows for {gene_id}-T1"
    return exons


def accepted_strains(path: str) -> set[str]:
    return set(l.strip() for l in open(path) if l.strip())


def count_exon_indels(exons, accepted) -> dict[str, object]:
    """For each gene: total indel records overlapping its CDS exons, and among those
    how many segregate (>=1 alt call) within the 213 accepted panel."""
    total_overlap = 0
    segregating = 0
    details = []
    for (scf, s, e) in exons:
        region = f"{scf}:{s}-{e}"
        r = subprocess.run([BCFTOOLS, "query", "-r", region,
                            "-f", "%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n", INDEL_VCF],
                           capture_output=True, text=True)
        for line in r.stdout.splitlines():
            parts = line.split("\t")
            pos, ref, alt = parts[0], parts[1], parts[2]
            total_overlap += 1
            seg = False
            for tok in parts[3:]:
                if "=" not in tok:
                    continue
                st, gt = tok.split("=")
                if st not in accepted:
                    continue
                allele = gt.replace("|", "/").split("/")[0]
                if allele not in (".", "0"):
                    seg = True
                    break
            if seg:
                segregating += 1
                details.append((pos, ref, alt))
    return dict(n_indels_in_cds=total_overlap, n_segregating_in_panel=segregating,
                segregating_sites=details)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-csv", required=True)
    args = ap.parse_args()

    accepted = accepted_strains(ACCEPTED)
    rows = []
    for gene in sorted(GENES):
        exons = load_cds_exons(GFF, gene)
        res = count_exon_indels(exons, accepted)
        rows.append(dict(gene_id=gene, n_exons=len(exons),
                         n_indel_records_in_cds=res["n_indels_in_cds"],
                         n_segregating_in_panel=res["n_segregating_in_panel"],
                         cds_indels_segregate=(res["n_segregating_in_panel"] > 0),
                         snp_only_premise_holds=(res["n_segregating_in_panel"] == 0)))
        print(f"{gene}: {res['n_indels_in_cds']} indel records in CDS exons, "
              f"{res['n_segregating_in_panel']} segregating in 213 panel "
              f"-> premise {'HOLDS' if res['n_segregating_in_panel']==0 else 'VIOLATED'}")
        if res["segregating_sites"]:
            print(f"   sites: {[(p, r, a) for (p, r, a) in res['segregating_sites'][:8]]}"
                  f"{' ...' if len(res['segregating_sites'])>8 else ''}")

    import pandas as pd
    df = pd.DataFrame(rows)
    pd.DataFrame(rows).to_csv(args.out_csv, index=False)
    print(f"\nWrote {args.out_csv}")
    print("Genes where the SNP-only premise does NOT hold:",
          df.loc[df["snp_only_premise_holds"] == False, "gene_id"].tolist())


if __name__ == "__main__":
    main()
