#!/usr/bin/env python3
"""scaffold_16 investigation (GWAS.md section 19): is the scaffold's recurring
GWAS-hit pattern an assembly/mapping artifact, one non-recombining haploblock,
or elevated background population differentiation?

Three checks, each against data already generated elsewhere in this project
(no new sequencing/genotyping):

1. Per-scaffold QC from the raw VCF (INFO/DP, INFO/AF) -- depth, SNP density,
   rare-variant fraction -- compared across all 23 nuclear scaffolds.
2. Pairwise LD (plink --r2) between the specific positions that recur as top
   hits for different traits on scaffold_16, on the full unpruned 213-strain
   marker set (the LD-PRUNED kinship set is unsuitable for this -- pruning
   deliberately removes correlated markers, so low r2 there is circular).
3. Per-scaffold mean Hudson Fst from the existing genome-wide pixy output.

Usage:
  pixi run python3 analysis/gwas/scripts/qc_scaffold16_investigation.py \
      --vcf data/raw/genotypes/RmucY2510_v2/RmucY2510_v2.All.SNP.combined_selected.vcf.gz \
      --bfile analysis/gwas/results/gwas/grm_conditioning/rebuilt_kinship/gwas \
      --pixy-fst analysis/ideas/2026-08-15-color-phenotype-space/results/gwas/pixy/genome_fst.txt \
      --bcftools <path> --plink <path> \
      --out-dir analysis/gwas/results/gwas/scaffold16_investigation
"""
from __future__ import annotations

import argparse
import pathlib
import subprocess

import pandas as pd

FOCAL_SNPS_BP = [121473, 122361, 416360, 418561, 455499, 492282, 563722]
RARE_AF_THRESHOLD = 0.05  # matches this project's rare-variant framing elsewhere (D-24/D-25)


def run_dp_af_scan(bcftools: str, vcf: str, out_dir: pathlib.Path) -> pathlib.Path:
    out = out_dir / "site_dp_af.tsv"
    with open(out, "w") as f:
        subprocess.run([bcftools, "query", "-f", "%CHROM\t%INFO/DP\t%INFO/AF\n", vcf],
                        check=True, stdout=f)
    return out


def per_scaffold_qc(site_dp_af_path: pathlib.Path, contig_lengths: dict[str, int]) -> pd.DataFrame:
    df = pd.read_csv(site_dp_af_path, sep="\t", header=None, names=["chrom", "dp", "af"])
    # bcftools query emits "." for sites where INFO/DP or INFO/AF is genuinely
    # absent (multiallelic AF edge cases); coercing to NaN is correct here --
    # rows still count toward n_sites (groupby key always present) and are
    # excluded only from the mean/dp,af aggregates themselves.
    df["af"] = pd.to_numeric(df["af"], errors="coerce")  # ANALYSIS_OK[missingness]: "." -> NaN, see comment above; still counted in n_sites
    df["dp"] = pd.to_numeric(df["dp"], errors="coerce")  # ANALYSIS_OK[missingness]: "." -> NaN, see comment above; still counted in n_sites
    g = df.groupby("chrom").agg(
        n_sites=("dp", "size"), mean_dp=("dp", "mean"),
        frac_af_lt_0p05=("af", lambda s: (s < RARE_AF_THRESHOLD).mean()),
    ).reset_index()
    g["length"] = g["chrom"].map(contig_lengths)
    g["sites_per_kb"] = g["n_sites"] / (g["length"] / 1000)
    return g.sort_values("sites_per_kb", ascending=False)


def contig_lengths_from_vcf_header(bcftools: str, vcf: str) -> dict[str, int]:
    out = subprocess.run([bcftools, "view", "-h", vcf], check=True, capture_output=True, text=True)
    lengths = {}
    for line in out.stdout.splitlines():
        if line.startswith("##contig=<ID="):
            parts = dict(kv.split("=", 1) for kv in line[len("##contig=<"):-1].split(",") if "=" in kv)
            lengths[parts["ID"]] = int(parts["length"])
    # ANALYSIS_OK[runtime-assert]: developer tripwire on the VCF header's
    # expected contig-line format, checked once at parse time.
    assert lengths, f"no ##contig lines parsed from {vcf} header"
    return lengths


def focal_snp_ld(plink: str, bfile: str, out_dir: pathlib.Path) -> pd.DataFrame:
    # ANALYSIS_OK[file-selection]: fixed plink1 output-suffix convention
    # (--bfile <prefix> always has exactly one <prefix>.bim), not a glob.
    bim = pd.read_csv(f"{bfile}.bim", sep=r"\s+", header=None,
                       names=["chrom", "snp_id", "cm", "bp", "a1", "a2"])
    focal = bim[(bim.chrom == 16) & (bim.bp.isin(FOCAL_SNPS_BP))]
    ids_path = out_dir / "focal_snp_ids.txt"
    ids_path.write_text("\n".join(focal.snp_id) + "\n")
    ld_out = out_dir / "focal_ld"
    subprocess.run([plink, "--bfile", bfile, "--allow-extra-chr", "--allow-no-sex",
                     "--extract", str(ids_path), "--r2", "--ld-window-r2", "0",
                     "--ld-window-kb", "700", "--ld-window", "999999",
                     "--out", str(ld_out)], check=True, capture_output=True, text=True)
    # ANALYSIS_OK[file-selection]: fixed plink1 --r2 --out output suffix, not a glob.
    return pd.read_csv(f"{ld_out}.ld", sep=r"\s+")


def pixy_fst_per_scaffold(pixy_fst_path: str) -> pd.DataFrame:
    df = pd.read_csv(pixy_fst_path, sep="\t")
    return (df.groupby("chromosome")["avg_hudson_fst"]
              .agg(mean_fst="mean", median_fst="median", n_windows="count")
              .sort_values("mean_fst", ascending=False).reset_index())


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--bfile", required=True, help="full unpruned 213-strain plink prefix")
    ap.add_argument("--pixy-fst", required=True)
    ap.add_argument("--bcftools", default="bcftools")
    ap.add_argument("--plink", default="plink")
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    out_dir = pathlib.Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("[scaffold16] 1/3 per-scaffold VCF QC (depth, density, rare-variant fraction) ...")
    lengths = contig_lengths_from_vcf_header(args.bcftools, args.vcf)
    dp_af_path = run_dp_af_scan(args.bcftools, args.vcf, out_dir)
    qc = per_scaffold_qc(dp_af_path, lengths)
    qc.to_csv(out_dir / "per_scaffold_qc.csv", index=False)  # ANALYSIS_OK[file-selection]: fixed output filename, not a glob
    print(qc.round(3).to_string(index=False))

    print("\n[scaffold16] 2/3 pairwise LD between recurring scaffold_16 hit positions ...")
    ld = focal_snp_ld(args.plink, args.bfile, out_dir)
    ld.to_csv(out_dir / "focal_snp_ld.csv", index=False)  # ANALYSIS_OK[file-selection]: fixed output filename, not a glob
    print(ld.to_string(index=False))

    print("\n[scaffold16] 3/3 per-scaffold mean Hudson Fst (existing pixy output) ...")
    fst = pixy_fst_per_scaffold(args.pixy_fst)
    fst.to_csv(out_dir / "per_scaffold_fst.csv", index=False)  # ANALYSIS_OK[file-selection]: fixed output filename, not a glob
    print(fst.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
