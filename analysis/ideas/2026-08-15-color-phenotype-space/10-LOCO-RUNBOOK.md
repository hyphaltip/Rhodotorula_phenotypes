# LOCO Sensitivity Runbook — Confirm Tier-A GWAS Sensitivity (D-10)

**Status**: ready to run (planned follow-up; apply after Tier A/B/C)
**See also**: `.living/decisions.md` D-10, `todo/loco-followup.md`, `.living/learnings.md` L-19

---

## Purpose

Tier-A next-gen GWAS scans (12 traits × all-201 + culled-173) use ONE full kinship GRM
(`gemma -lmm 4 -k output/kins.cXX.txt`). In an LMM, a SNP's own scaffold's genetic similarity
can be over-absorbed by the kinship term ("proximal contamination"), biasing true associations
toward under-detection. **LOCO** (leave-one-chromosome-out) rebuilds the GRM excluding the
scaffold under test so each SNP is tested against a kinship built only from the OTHER scaffolds.

The user explicitly wants LOCO applied to **confirm we had appropriate sensitivity** before
finalizing conclusions — this is a committed sensitivity check, not an optimization of the
initial scan.

---

## Prerequisites (all exist / reusable)

| Item | Path | Notes |
|------|------|-------|
| GEMMA 0.98.3 | `/bigdata/stajichlab/shared/projects/Population_Genomics/Rhodotorula_mucilaginosa_NRRLY2510/GWAS/vc2gwas_env/bin/gemma` | no `-loco` flag → per-scaffold approach (L-19) |
| plink2 | `/opt/linux/centos/8.x/x86_64/pkgs/plink2/2.0.0-a.7.3/bin/intel_avx2/plink2` (path varies by node) | for per-scaffold variant subsetting |
| Trait bfiles (full) | `/scratch/jstajich/27384933/gwas/work/ngwas_<t>.{bed,bim,fam}` | trait in .fam col-6; .bim integer chr 1..23 |
| Culled trait bfiles | `/scratch/jstajich/27384933/gwas/work/ngwasc_<t>.{bed,bim,fam}` | 173 strains; trait in col-6 |
| Pruned kinship variants (full) | `gwas.pruned.{bed,bim,fam}` | 20,769 markers, integer chr codes |
| Pruned kinship variants (culled) | `gwasc.pruned.{bed,bim,fam}` | 173 strains |
| LOCO wrapper | `/scratch/jstajich/27384933/gwas/work/run_loco.sh` | VALIDATED on scaffold_1/chroma |

CHR inventory in pruned set: scaffold_1..23 → 1..23 (counts per chr range 5–2,764; scaffold_6=74,
11=47, 18=34, 20=5).

---

## How LOCO works here (GEMMA 0.98.3 — no -loco flag)

For each scaffold `c`:
1. Take the **pruned** variant set and **drop scaffold `c`** → build that leave-out GRM
   (20,769 minus chr-`c` markers).
2. Take the **full trait bfile** and extract **only scaffold `c`'s polymorphic SNPs** for testing.
3. Run `gemma -lmm 4 -k <leaveout_GRM_c>` on scaffold-`c` SNPs.
   Result: each SNP on `c` is tested against a kinship built only from other scaffolds.

Repeat for all scaffolds; concatenate per-scaffold `.assoc.txt` into a genome-wide LOCO scan.

---

## Run (single warp per trait, serial inside)

```bash
cd /scratch/jstajich/27384933/gwas/work
# all-201, e.g. chroma
bash run_loco.sh chroma ngwas_chroma gwas.pruned ngwas 1        # 'ngwas' = outroot prefix
# culled-173, e.g. chroma
bash run_loco.sh chroma ngwasc_chroma gwasc.pruned ngwasc 1
```

`run_loco.sh` args: `<trait> <bfile-root> <kinship-root> <outroot> [job] [scaffold]`.

Outputs: `output/loco_<outroot>_<trait>_<chr>.assoc.txt` per scaffold.

### Optional SLURM array (one task per scaffold) — recommended for all-23 across traits
```bash
#SBATCH --array=1-23
scaff=$(cut -f1 "gwas.pruned.bim" | sort -n | uniq | sed -n "${SLURM_ARRAY_TASK_ID}p")
bash run_loco.sh chroma ngwas_chroma gwas.pruned ngwas "$SLURM_ARRAY_TASK_ID" "$scaff"
```
Submit from a SHARED workdir (`sbatch -D /bigdata/...`), not node-local scratch (L-18).

---

## Recommended procedure

1. Wait for Tier A (24 scans) + Tier B (SKAT) + Tier C (BSLMM) to finish.
2. **Pick the top 1–3 traits** by strongest Tier-A signal (and Tier-B SKAT hits).
3. Run LOCO on those traits, all-201 **and** culled-173 (6 runs max).
4. For each trait: **merge** per-scaffold LOCO assoc → genome-wide file; recount lambda, top
   p-values, and independent loci after clump.
5. **Compare** LOCO vs full-kinship Tier-A:
   - lambda / PVE
   - top-hit identity and rank
   - whether any chromosome-proximal top hit gains or loses significance
6. Document the comparison in `GWAS_REPORT.md` Stage 6 (see acceptance criteria below).

---

## Interpretation & acceptance criteria

- **Confirmed sensitivity**: top hits and lambda stable between full-kinship and LOCO → Tier-A
  conclusions hold; proximal contamination did not mask signal.
- **Sensitivity gap found**: any trait where a scaffold-proximal signal moves substantially
  (e.g. a top hit appears/gains after LOCO) → report that full-kinship under-detected it; report
  the LOCO result for that trait as the more sensitive estimate.
- Either way, write a one-line conclusion in the report: whether the full-kinship Tier-A
  under-detected chromosome-proximal loci.

---

## Gotchas (learned)

- GEMMA 0.98.3 has **no** `-loco` flag (L-19) — the per-scaffold GRM loop is the only way.
- Leave-out GRM must be built from **all strains** (same order) so `-k` rows match the `.fam`;
  the wrapper copies the trait bfile `.fam` onto the per-scaffold testset.
- GEMMA kinship `.fam` col-6 must be non-`-9`; the wrapper zeroes it for `-gk 1`.
- Scaffolds with very few pruned markers (scaffold_6=74, 11=47, 18=34, 20=5) still produce valid
  GRMs but with little leave-out differential — note this in per-scaffold logs.
- plink2 path is **node-dependent** (`/opt/linux/centos/...` on some nodes vs `/opt/linux/rocky/...`).