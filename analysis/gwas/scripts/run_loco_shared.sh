#!/bin/bash
# Shared-storage LOCO (leave-one-chromosome-out) wrapper for GEMMA 0.98.3.
# Variant of run_loco.sh hardened for SLURM batch: WORK = shared /bigdata dir, so
# batch nodes (which cannot see node-local scratch) can reach bfiles. (L-18)
#
# Usage:  bash run_loco_shared.sh <trait> <bfile-root> <kinship-root> <outroot> [job] [scaff]
#   WORK    env or arg7: shared work dir containing bfiles/{bed,bim,fam} + pruned + plink_arch.sh
#   <bfile-root> bfile basename WITHOUT ext (trait in .fam col-6, .bim integer chr 1..23)
#   <kinship-root> pruned bfile root (integer chr codes) used to build leave-out GRMs
#   <outroot> output prefix; assoc written to $WORK/output/loco_<outroot>_<trait>_<chr>.assoc.txt
#   <scaff>  restrict to one scaffold (SLURM array tasks) or empty for all
set -euo pipefail

WORK="${7:?WORK dir required as arg7 (must contain plink_arch.sh + bfiles)}"
cd "$WORK"
GEMMA=/bigdata/stajichlab/shared/projects/Population_Genomics/Rhodotorula_mucilaginosa_NRRLY2510/GWAS/vc2gwas_env/bin/gemma

TRAIT="${1:?need trait}"
BFILE="${2:?need bfile-root}"
KINROOT="${3:?need kinship-root (pruned)}"
OUTROOT="${4:?need outroot}"
JOBID="${5:-1}"
SCAFF="${6:-}"

mkdir -p "$WORK/output"
echo "### LOCO trait=$TRAIT bfile=$BFILE kinroot=$KINROOT [scaff=$SCAFF] $(date) ###"

CHRS=$(cut -f1 "$KINROOT.bim" | sort -n | uniq)
[ -n "$SCAFF" ] && CHRS="$SCAFF"

for c in $CHRS; do
  echo "=== scaffold $c ==="
  awk -v c=$c '$1 != c {print $2}' "$KINROOT.bim" > /tmp/loco_varids_${c}_${JOBID}.txt
  [ -n "${PLINK:-}" ] || source "$WORK/plink_arch.sh"
  $PLINK --bfile "$KINROOT" --extract /tmp/loco_varids_${c}_${JOBID}.txt --make-bed --out /tmp/loco_${c}_grmset_${JOBID} 2>&1 | tail -1
  awk '{print $1,$2,$3,$4,$5,0}' /tmp/loco_${c}_grmset_${JOBID}.fam > /tmp/loco_${c}_grmset_${JOBID}.fam0
  mv /tmp/loco_${c}_grmset_${JOBID}.fam0 /tmp/loco_${c}_grmset_${JOBID}.fam
  $GEMMA -bfile /tmp/loco_${c}_grmset_${JOBID} -gk 1 -o loco_${c}_kin_${JOBID} 2>&1 | tr '\r' '\n' | tail -1

  awk -v c=$c '$1==c{print $2}' "$BFILE.bim" > /tmp/loco_test_${c}_${JOBID}.txt
  [ -s /tmp/loco_test_${c}_${JOBID}.txt ] || { echo "no variants on scaffold $c in bfile; skip"; continue; }
  $PLINK --bfile "$BFILE" --extract /tmp/loco_test_${c}_${JOBID}.txt --make-bed --out /tmp/loco_${c}_testset_${JOBID} 2>&1 | tail -1
  cp "$BFILE.fam" /tmp/loco_${c}_testset_${JOBID}.fam
  $GEMMA -bfile /tmp/loco_${c}_testset_${JOBID} -k "$WORK/output/loco_${c}_kin_${JOBID}.cXX.txt" -lmm 4 \
         -o loco_${OUTROOT}_${TRAIT}_${c} 2>&1 | tr '\r' '\n' | tail -1

  rm -f /tmp/loco_${c}_grmset_${JOBID}.{bed,bim,fam,nosex,log} /tmp/loco_${c}_testset_${JOBID}.{bed,bim,fam,nosex,log} \
        /tmp/loco_varids_${c}_${JOBID}.txt /tmp/loco_test_${c}_${JOBID}.txt /tmp/loco_${c}_grmset_${JOBID}*
  rm -f "$WORK/output/loco_${c}_kin_${JOBID}."*
done
echo "LOCO_DONE trait=$TRAIT outroot=$OUTROOT (scaffolds: $CHRS)"
