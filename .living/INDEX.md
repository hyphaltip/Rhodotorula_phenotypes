# .living/ Index
Last audit: 2026-08-16

| File | Entries | Last updated | Key topics |
|------|---------|--------------|------------|
| conventions.md | 0 sections | 2026-08-15 | — |
| decisions.md | 9 entries | 2026-08-16 | D-9 — Use kinship-only GEMMA LMM as primary GWAS; park PC/pop-covariate models (singular GRM from near-clone strains), D-4 — Make strain metadata authoritative over CSV + add `--strain-only` sync to the DB import, D-5 — Serialize repeated analyses through a shared extract (db_extract) instead of re-querying DuckDB per idea script, D-6 — Pin numpy to 2.4.x in pixi and disable user-site so the conda env is authoritative, D-7 — Use Mantel permutation as the primary phylogenetic-signal statistic; keep K/lambda as caveated references, D-3 — Stratify plate-effect variance partition by Copper concentration |
| learnings.md | 18 entries | 2026-08-16 | L-18 — UCR HPCC batch jobs don't inherit login shell: bgzip/tabix/pixi unavailable; submit from shared cwd (node-local scratch kills job, signal 53), L-17 — Structure covariates (PCs/pop dummies) collapse GEMMA LMM when kinship singular from near-clone redundancy, L-16 — GEMMA reads phenotype from fam col-6, silently ignores `-p` w/ -bfile; needs integer chrom codes, L-1 — mycelium scripts require Python >= 3.10 |
| log/ | 6 sessions | 2026-08-16 | rhodotorula-phenotypes (6) |
| findings/ | 2 findings across 2 topics | 2026-08-16 | phenotype-color-space, copper-tolerance |

## Local skills
See `.living/skills/` for project-specific skill packs.
