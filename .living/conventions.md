# Repo-Specific Conventions

Overrides to mycelium defaults or convention pack conventions.

<!-- Document any project-specific convention overrides here. -->

## GEMMA/GWAS pipeline conventions (this repo's `analysis/gwas/` and prior `analysis/ideas/2026-08-15-color-phenotype-space/`)

Distilled from 11+ recurring gotchas across two GWAS analyses (L-16 through L-26). Follow these whenever writing or rerunning GEMMA-based GWAS scripts in this repo:

1. **`.bim` chromosome codes must be integers** (GEMMA requires it) — remap `scaffold_N` -> `N` via `Chrom_Mapping.tab` before any `-gk`/`-lmm` call, or GEMMA silently reports "0 analyzed individuals".
2. **GEMMA reads the phenotype from `.fam` column 6, not from `-p`, whenever `-bfile` and a kinship/`-p` are both given.** Bake the trait into a per-trait `.fam` copy (missing -> `NA`) rather than relying on `-p`.
3. **Kinship (GRM) and the association scan use DIFFERENT SNP sets, by design**: LD-pruned SNPs (`--indep-pairwise 50 5 0.2`) build the kinship matrix; the FULL QC'd unpruned SNP set is scanned for association. Persist BOTH bfiles under unambiguous, differently-named files in the same rebuild step — never let the unpruned set live only on ephemeral `$SCRATCH` (L-26: this exact mistake silently regressed a rebuild to pruned-only association, testing <10% of the genome).
4. **Never use `$(dirname "${BASH_SOURCE[0]}")` in a script destined for a SLURM job** — resolves to nothing/wrong path on batch nodes. Use `$PWD` with an explicit "run from repo root" assertion, or pass paths as arguments/env vars.
5. **This cluster's plink2 build (`vc2gwas_env/bin/plink2`) has no `--distance` flag.** Use `--make-king-table cols=id,nsnp,ibs0` for pairwise IBS/kinship needs instead (L-25).
6. **BH-FDR, not permutation, for significance** on this near-clonal panel — a manual maxT phenotype-permutation null is uncalibrated here (fixed-kinship + permuted-phenotype variance outliers collapse the tail; D-11/L-21). Report a Meff-proxy (LD-pruned SNP count) + Bonferroni-at-Meff alongside BH-FDR as sensitivity context, not a replacement (D-16).
7. **Eigendecompositions of real-data LD/kinship matrices need a shift-retry + scipy fallback** (`_safe_eigvals` pattern) — near-singular matrices routinely fail plain `np.linalg.eigvalsh` (L-22).
8. **Any algorithm described only in a PROGRESS.md prose table (not saved as a script) must be reconstructed AND validated against the original's known output before trusting it on new data** — e.g. the near-clone IBS0-culling algorithm (D-15/D-16): reconstructed, then checked for >90% membership overlap against the prior run's known 173-strain result before applying to a new panel.

Source: L-16, L-17, L-18, L-19, L-20, L-21, L-22, L-23, L-24, L-25, L-26; D-9, D-11, D-15, D-16.

## HPC module-load convention (this cluster, applies to any script/interactive shell)

**Never pipe a `module load ...` invocation through anything** (`| tail`, `| grep`, `2>&1 | ...`). Piping forces a subshell, and `module` is a shell function that mutates the current shell's environment via `eval` — a subshell's mutation vanishes the instant the pipe exits, so the load silently "succeeds" (no error surfaces) while every env var/PATH change it should have made never lands in the calling shell. Always run `module load foo` as a bare statement, in the same shell as the commands depending on it; verify with a separate `env | grep -i foo` if needed, never chained onto the load itself with a pipe.

Recurred 3x in one session (hmmer, kofamscan, snpEff) before being crystallized here — see L-30.

Source: L-30.

## Convention-pack feedback: robust-analysis on the candidate-gene alignment pilot (2026-08-26)

The `robust-analysis` core convention's "fail loudly on unexpected data / assert shapes-types-ranges" practice paid off directly in `analysis/candidate_gene_alignment/scripts/extract_gene_sequences.py`: the hard assertion that every VCF REF allele matches the genome base at build time, and that all 213 strains produce equal-length CDS, are exactly what let the NRRL_Y-2510 reference-strain sanity check (D-21) be trusted as a real correctness check rather than a coincidence — if either assertion had been silently skipped, a coordinate or REF/ALT-swap bug could have passed unnoticed. No gaps found in the convention for this task; it directly caught the class of bug it's designed to catch, even though no bug was actually present this time.

Source: D-21.
