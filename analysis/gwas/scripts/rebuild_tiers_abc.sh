#!/usr/bin/env bash
# Rebuild kinship + GEMMA Tier A/B/C on a new strain set. Only run when
# state_diff_report.json says overall_action == "rebuild" (Task 5).
# Reconstructed from analysis/ideas/2026-08-15-color-phenotype-space/PROGRESS.md
# sections 1-2 (S0-S4) -- the original commands ran ad hoc on $SCRATCH and were
# never saved as a reusable script until this port.
#
# NOTE: intentionally does not resolve its own path via BASH_SOURCE (breaks on
# SLURM work directories) -- run this via `sbatch --chdir=<repo_root> --wrap=...`
# or `cd <repo_root> && sbatch --wrap="bash analysis/gwas/scripts/rebuild_tiers_abc.sh"`.
set -euo pipefail
module load bcftools

REPO_ROOT="$PWD"
if [ ! -f "$REPO_ROOT/pixi.toml" ] || [ ! -d "$REPO_ROOT/analysis/gwas" ]; then
  echo "ERROR: expected to be run from the repo root (pixi.toml + analysis/gwas/ not found in $REPO_ROOT)." >&2
  echo "Re-run as: cd <repo_root> && sbatch --wrap=\"bash analysis/gwas/scripts/rebuild_tiers_abc.sh\"" >&2
  exit 1
fi
cd "$REPO_ROOT"

: "${SCRATCH:?SCRATCH env var not set -- run this inside a SLURM job}"
WORK="$SCRATCH/gwas_rebuild"
mkdir -p "$WORK"

VCF="data/raw/genotypes/RmucY2510_v2/RmucY2510_v2.All.SNP.combined_selected.vcf.gz"
STRAIN_LIST="analysis/gwas/results/strain_reconciliation/accepted_vcf_ids.txt"
TOOLCHAIN="/bigdata/stajichlab/shared/projects/Population_Genomics/Rhodotorula_mucilaginosa_NRRLY2510/GWAS/vc2gwas_env/bin"

# Emit the accepted strain list (post-reconciliation, post-ploidy-exclusion) as a
# plain VCF-sample-ID list for bcftools -S.
pixi run python3 -c "
import csv
with open('analysis/gwas/results/strain_reconciliation/strain_match_table.reviewed.csv') as f:
    ids = [r['vcf_sample_id'] for r in csv.DictReader(f)
           if (r.get('tier') in ('exact','normalized') or r.get('decision') == 'accept') and r['vcf_sample_id']]
with open('$STRAIN_LIST', 'w') as out:
    out.write('\n'.join(sorted(set(ids))) + '\n')
print(f'{len(set(ids))} accepted strain IDs written to $STRAIN_LIST')
"

# S0: subset
bcftools view -S "$STRAIN_LIST" "$VCF" -Oz -o "$WORK/sub.vcf.gz"
bcftools index -t "$WORK/sub.vcf.gz"

# S2: biallelic
bcftools view -m2 -M2 "$WORK/sub.vcf.gz" -Oz -o "$WORK/sub.biallelic.vcf.gz"

# S3: plink2 QC + LD-pruned set for kinship (mirrors PROGRESS.md S3/S3b thresholds)
# The source VCF has "." for every variant ID (bcftools view confirms this) -- plink2's
# --indep-pairwise requires unique IDs, so assign chrom:pos:ref:alt IDs at import time
# (undocumented in PROGRESS.md's ad hoc run; this port makes it explicit and reusable).
"$TOOLCHAIN/plink2" --vcf "$WORK/sub.biallelic.vcf.gz" --set-all-var-ids '@:#:$r:$a' \
  --new-id-max-allele-len 200 --geno 0.1 --mind 0.1 --mac 5 \
  --make-bed --out "$WORK/gwas" --allow-extra-chr
"$TOOLCHAIN/plink2" --bfile "$WORK/gwas" --indep-pairwise 50 5 0.2 \
  --out "$WORK/gwas.pruned" --allow-extra-chr
"$TOOLCHAIN/plink2" --bfile "$WORK/gwas" --extract "$WORK/gwas.pruned.prune.in" \
  --make-bed --out "$WORK/gwas.pruned" --allow-extra-chr

# GEMMA needs integer chromosome codes (learning L-16) -- rename per the source
# project's chrom map before running -gk/-lmm.
CHROM_MAP="/bigdata/stajichlab/shared/projects/Population_Genomics/Rhodotorula_mucilaginosa_NRRLY2510/genome/Chrom_Mapping.tab"
pixi run python3 -c "
import pandas as pd
m = pd.read_csv('$CHROM_MAP', sep='\t', header=None, names=['scaffold','int_chr'], skiprows=1)
mapping = dict(zip(m['scaffold'], m['int_chr']))
bim = pd.read_csv('$WORK/gwas.pruned.bim', sep='\t', header=None)
bim[0] = bim[0].map(mapping)
assert bim[0].notna().all(), 'unmapped scaffolds in .bim after chrom-code remap'
bim.to_csv('$WORK/gwas.pruned.bim', sep='\t', header=False, index=False)
"

# Kinship (needs -p even for -gk 1, per L-16 gotcha)
"$TOOLCHAIN/gemma" -bfile "$WORK/gwas.pruned" -p "$WORK/gwas.pruned.fam" -gk 1 \
  -o kins -outdir "$WORK/output"

# GRM conditioning diagnostic -- MANDATORY on any strain-state change (spec S4)
mkdir -p analysis/gwas/results/gwas/grm_conditioning
pixi run python3 analysis/gwas/scripts/check_grm_conditioning.py \
  --kinship "$WORK/output/kins.cXX.txt" \
  --out analysis/gwas/results/gwas/grm_conditioning/grm_diagnostic.json

# Persist the pruned .fam/.bed/.bim and kinship matrix back to shared storage --
# $SCRATCH is node-local and will not survive past this job.
mkdir -p analysis/gwas/results/gwas/grm_conditioning/rebuilt_kinship
cp "$WORK/gwas.pruned.fam" "$WORK/gwas.pruned.bed" "$WORK/gwas.pruned.bim" \
   "$WORK/output/kins.cXX.txt" \
   analysis/gwas/results/gwas/grm_conditioning/rebuilt_kinship/

echo "Kinship + GRM diagnostic complete. Rebuilt kinship + pruned genotype set copied to"
echo "analysis/gwas/results/gwas/grm_conditioning/rebuilt_kinship/. Next: rerun"
echo "build_gwas_phenotypes.py + per-trait 'gemma -lmm 4 -k' (Tier A), then"
echo "tierb_set_tests.py (Tier B) and the BSLMM invocation (Tier C) against this"
echo "rebuilt kinship, following PROGRESS.md sections 6-9 for exact per-trait commands."
