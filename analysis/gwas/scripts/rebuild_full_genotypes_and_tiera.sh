#!/usr/bin/env bash
# Rebuild the FULL (unpruned, QC'd) genotype sets for both the all-213 ("gwas")
# and near-clone-culled 182-strain ("gwasc") panels, persist them to shared
# storage (Task 7's rebuild_tiers_abc.sh only persisted the LD-pruned kinship
# set -- Tier A association testing needs the full QC'd marker set, per the
# original run: 404,706 SNPs for association vs. 20,769 for kinship. This was
# a scope gap in the first Tier A rebuild -- fixed here), build per-panel
# kinship (already have gwas's from Task 7; builds gwasc's here), and rerun
# Tier A (12 traits) on BOTH panels against the FULL marker set.
#
# Run from the repo root, inside a SLURM job (needs $SCRATCH):
#   sbatch --wrap="bash analysis/gwas/scripts/rebuild_full_genotypes_and_tiera.sh"
set -euo pipefail
module load bcftools

REPO_ROOT="$PWD"
[ -f "$REPO_ROOT/pixi.toml" ] || { echo "ERROR: run from repo root." >&2; exit 1; }

: "${SCRATCH:?SCRATCH env var not set -- run inside a SLURM job}"
WORK="$SCRATCH/gwas_full_rebuild"
mkdir -p "$WORK"

VCF="data/raw/genotypes/RmucY2510_v2/RmucY2510_v2.All.SNP.combined_selected.vcf.gz"
STRAIN_LIST="analysis/gwas/results/strain_reconciliation/accepted_vcf_ids.txt"
CULLED_LIST="analysis/gwas/results/gwas/near_clone_culling/culled_keep.txt"
TOOLCHAIN="/bigdata/stajichlab/shared/projects/Population_Genomics/Rhodotorula_mucilaginosa_NRRLY2510/GWAS/vc2gwas_env/bin"
CHROM_MAP="/bigdata/stajichlab/shared/projects/Population_Genomics/Rhodotorula_mucilaginosa_NRRLY2510/genome/Chrom_Mapping.tab"
PHENO_CSV="analysis/gwas/results/gwas_next_phenotypes_fam_order.csv"

for f in "$VCF" "$STRAIN_LIST" "$CULLED_LIST" "$PHENO_CSV"; do
  [ -f "$f" ] || { echo "ERROR: required input missing: $f" >&2; exit 1; }
done

OUTDIR_BASE="analysis/gwas/results/gwas"
mkdir -p "$OUTDIR_BASE/grm_conditioning/rebuilt_kinship" "$OUTDIR_BASE/tierA_summary/gemma_output"

remap_chrom() {
  # $1 = .bim path to remap in place (scaffold_N -> integer N per Chrom_Mapping.tab)
  pixi run python3 -c "
import pandas as pd
m = pd.read_csv('$CHROM_MAP', sep='\t', header=None, names=['scaffold','int_chr'], skiprows=1)
mapping = dict(zip(m['scaffold'], m['int_chr']))
bim = pd.read_csv('$1', sep='\t', header=None)
bim[0] = bim[0].map(mapping)
assert bim[0].notna().all(), 'unmapped scaffolds in $1 after chrom-code remap'
bim.to_csv('$1', sep='\t', header=False, index=False)
"
}

# ---------------------------------------------------------------------------
# Panel: gwas (all 213, unpruned QC'd set)
# ---------------------------------------------------------------------------
echo "=== gwas panel: subset + QC (unpruned) ==="
bcftools view -S "$STRAIN_LIST" "$VCF" -Oz -o "$WORK/sub.vcf.gz"
bcftools index -t "$WORK/sub.vcf.gz"
bcftools view -m2 -M2 "$WORK/sub.vcf.gz" -Oz -o "$WORK/sub.biallelic.vcf.gz"

"$TOOLCHAIN/plink2" --vcf "$WORK/sub.biallelic.vcf.gz" --set-all-var-ids '@:#:$r:$a' \
  --new-id-max-allele-len 200 --geno 0.1 --mind 0.1 --mac 5 \
  --make-bed --out "$WORK/gwas" --allow-extra-chr
remap_chrom "$WORK/gwas.bim"
cp "$WORK/gwas.bed" "$WORK/gwas.bim" "$WORK/gwas.fam" "$OUTDIR_BASE/grm_conditioning/rebuilt_kinship/"
echo "gwas: $(wc -l < "$WORK/gwas.bim") variants x $(wc -l < "$WORK/gwas.fam") strains"

# ---------------------------------------------------------------------------
# Panel: gwasc (182 near-clone-culled, unpruned QC'd set + pruned kinship)
# ---------------------------------------------------------------------------
echo "=== gwasc panel: culled-182 subset + QC (unpruned) + kinship ==="
awk '{print "0",$1}' "$CULLED_LIST" > "$WORK/culled.keep.txt"
"$TOOLCHAIN/plink2" --bfile "$WORK/gwas" --keep "$WORK/culled.keep.txt" \
  --make-bed --out "$WORK/gwasc" --allow-extra-chr
cp "$WORK/gwasc.bed" "$WORK/gwasc.bim" "$WORK/gwasc.fam" "$OUTDIR_BASE/grm_conditioning/rebuilt_kinship/"
echo "gwasc: $(wc -l < "$WORK/gwasc.bim") variants x $(wc -l < "$WORK/gwasc.fam") strains"

# LD-prune gwasc for its own kinship (mirrors Task 7's gwas.pruned derivation)
"$TOOLCHAIN/plink2" --bfile "$WORK/gwasc" --indep-pairwise 50 5 0.2 \
  --out "$WORK/gwasc.pruned" --allow-extra-chr
"$TOOLCHAIN/plink2" --bfile "$WORK/gwasc" --extract "$WORK/gwasc.pruned.prune.in" \
  --make-bed --out "$WORK/gwasc.pruned" --allow-extra-chr
"$TOOLCHAIN/gemma" -bfile "$WORK/gwasc.pruned" -p "$WORK/gwasc.pruned.fam" -gk 1 \
  -o kinsc -outdir "$WORK/output"
mkdir -p "$OUTDIR_BASE/grm_conditioning/rebuilt_kinship_culled"
cp "$WORK/gwasc.pruned.bed" "$WORK/gwasc.pruned.bim" "$WORK/gwasc.pruned.fam" \
   "$WORK/output/kinsc.cXX.txt" "$OUTDIR_BASE/grm_conditioning/rebuilt_kinship_culled/"

# GRM conditioning diagnostic for the culled panel too (spec S4 -- report, don't skip)
pixi run python3 analysis/gwas/scripts/check_grm_conditioning.py \
  --kinship "$WORK/output/kinsc.cXX.txt" \
  --out "$OUTDIR_BASE/grm_conditioning/grm_diagnostic_culled.json"

# ---------------------------------------------------------------------------
# Tier A: 12 traits x {gwas full-213, gwasc culled-182}, full unpruned SNP set,
# association against each panel's OWN pruned-set kinship (matches original
# design: kinship from pruned set, association test on full QC'd set).
# ---------------------------------------------------------------------------
TRAITS="chroma sat bright clone_mean_area AUC_0 AUC_10 AUC_20 AUC_30 AUC_ratio_10 resilience_30 cu_dose_slope IC50_est"

run_tiera_panel() {
  local panel="$1" bfile="$2" kinship="$3"
  for trait in $TRAITS; do
    echo "=== Tier A [$panel]: $trait ==="
    TWORK="$WORK/tiera_${panel}_${trait}"
    mkdir -p "$TWORK"
    cp "$bfile.bed" "$TWORK/g.bed"
    cp "$bfile.bim" "$TWORK/g.bim"
    pixi run python3 -c "
import pandas as pd
fam = pd.read_csv('$bfile.fam', sep=r'\s+', header=None, names=['fid','iid','pid','mid','sex','pheno'])
pheno = pd.read_csv('$PHENO_CSV')
pheno_map = dict(zip(pheno['strain_code'], pheno['$trait']))
vals = fam['iid'].map(pheno_map)
n_missing = vals.isna().sum()
fam['pheno'] = vals.apply(lambda v: 'NA' if pd.isna(v) else v)
fam.to_csv('$TWORK/g.fam', sep=' ', header=False, index=False)
print(f'  [$panel] $trait: {len(fam) - n_missing}/{len(fam)} with phenotype, {n_missing} NA')
"
    "$TOOLCHAIN/gemma" -bfile "$TWORK/g" -k "$kinship" -lmm 4 \
      -o "${panel}_${trait}" -outdir "$OUTDIR_BASE/tierA_summary/gemma_output" \
      > "$OUTDIR_BASE/tierA_summary/gemma_output_${panel}_${trait}.log" 2>&1 \
      || { echo "ERROR: GEMMA failed for $panel/$trait" >&2; exit 1; }
    rm -rf "$TWORK"
  done
}

run_tiera_panel gwas "$OUTDIR_BASE/grm_conditioning/rebuilt_kinship/gwas" \
  "$OUTDIR_BASE/grm_conditioning/rebuilt_kinship/kins.cXX.txt"
run_tiera_panel gwasc "$OUTDIR_BASE/grm_conditioning/rebuilt_kinship_culled/gwasc" \
  "$OUTDIR_BASE/grm_conditioning/rebuilt_kinship_culled/kinsc.cXX.txt"

echo "Full-marker-set Tier A rebuild complete: 12 traits x 2 panels."
