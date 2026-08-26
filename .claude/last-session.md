# Last session: 2026-08-26 (005, continued) — Candidate-gene alignment: 8-gene extension, pilot indel-screen correction, and color-phenotype association

## Done

- **Extended the candidate-gene alignment pipeline to all 8 GWAS-locus genes** via the
  new `scripts/run_extend_to_gwas_genes.sh` (CDS/protein FASTA + indel screen + snpEff
  region + variant table + alignment PNG, all in `results/OM429_{001533,001415,001430,
  003729,002663,005034,001521,000065}/`).

- **Found and corrected a silent-0 bug in the pilot's indel screen.** `screen_indels.sh`
  had reported 0 indels for both carotenoid genes; it was actually counting an empty
  bcftools pipe (tool not on PATH, `2>/dev/null` + `|| true`). Real counts: 433 and 821
  records in those CDS±2kb; `check_gene_coding_indels.py` shows **segregating in-CDS
  indels in 10/10 target genes**. Fixed script to fail loudly (L-31).

- **Added the color-phenotype association step** with the same population-aware battery
  the GWAS port validated (within-pop re-test + meta + covariate + FDR):
  - `associate_variants_with_phenotype.py`: 954 SNP-coding tests → 39 FDR-sig /
    40 `likely_real`.
  - `associate_indels_with_phenotype.py`: 2046 in-CDS-indel tests → 32 FDR-sig.
  - Results in `results/candidate_{gene,indel}_phenotype_assoc_all10.csv`.

- **Headline finding**: `OM429_000065` (sat-locus gene) c.839C>T p.Ala280Val (the lead
  SNP 208569) **and** a +6 bp in-frame GAGCGG-repeat insertion at 208398 in perfect LD
  with it — both associate with lab_a/chroma/sat (meta_p ≈ 8.8e-11, FDR ≈ 1e-8, partial
  R² 0.086–0.10, replicated 3/3 testable populations). The insertion is invisible to
  SNP-only pipelines (L-33). Other replicated: OM429_001521 p.Gln245His, OM429_005034
  p.Tyr12His, OM429_002663 p.Asn4Lys, OM429_001533 p.Tyr427Phe.

- **OM429_003333 splice_donor is monomorphic (0/213 alt)** → the "check splice_donor vs
  color" TODO resolves to N.A.; no testable variation.

- **Carotenoid genes (OM429_003333/003336) show no replicated color association** on
  either SNP or indel surface — consistent with no co-localization with validated loci.

## Decisions
- D-22 — correct indel screen + extend to 8 GWAS-locus genes + population-aware
  association step (SNP and indel).

## Learnings
- L-31 silent-0 screen toolchain failure; L-32 background-shell toolchain fixes (abs
  PATH / LD_LIBRARY_PATH / snpEff dataDir); L-33 coding repeat-copy variant carries the
  sat-locus signal, invisible to SNP-only analysis.

## Next steps
- Functional/expression validation of the sat-locus +6 bp GAGCGG insertion vs the
  lead SNP (perfect LD; cannot split causal without denser data) — todo added.
- Interactive alignment viewer (deferred per user scoping).
- Confirmatory annotation (KofamScan/BLASTp) of the nominated coding variants.
