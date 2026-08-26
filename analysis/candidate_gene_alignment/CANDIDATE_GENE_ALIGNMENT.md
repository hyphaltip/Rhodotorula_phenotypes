# Candidate-gene DNA/protein sequence alignment across the 213-strain panel

**Design doc**: `docs/superpowers/specs/2026-08-26-candidate-gene-alignment-design.md`
(reviewed by an Opus general-review pass and a bioinformatics-expert/fable consult
before implementation; both reviews' findings are addressed below).

**Status**: complete for all 10 target genes — 2 confirmed carotenoid pathway genes +
8 GWAS-locus nearest genes. Per-gene sequence alignment + variant tables + static
alignment PNGs all produced, PLUS the deferred phenotype-association step
(SNP coding variants AND indel genotypes vs 6 color traits), reported below.

## What this is (and isn't)

For each candidate gene, this builds a per-strain spliced CDS + protein sequence
across all 213 panel strains, a positional variant table (snpEff consequence/codon/AA
change per site, joined with each strain's genotype, population, phenotype value, and
— for GWAS-locus genes — the locus's own lead-SNP genotype and each coding variant's
own Tier A marginal p-value), and a static color-block alignment image.

The pilot described this as "descriptive, not a second association test." **That caveat
is now lifted with the phenotype-association step** (`scripts/
associate_variants_with_phenotype.py` + `associate_indels_with_phenotype.py`): coding
variants *are* formally tested against the 6 color traits, but only with the
population-aware battery (within-population re-test + fixed-effect meta-analysis +
population-covariate regression + BH-FDR), exactly the battery the GWAS port's
disambiguation exercise validated — NOT a naive pooled regression, which is
population-confounded on this near-clonal, structured panel.

## Correction to the pilot (important — read first)

The pilot's per-gene indel screen (`screen_indels.sh`) reported **0 indels** for both
carotenoid genes (and this was treated as confirming the "SNP-only substitution premise"
holds). That was **wrong**: the script silently returned 0 because bcftools was not on
PATH in that session (piped `2>/dev/null` + `grep -vc` + `|| true` swallowed the toolchain
failure). Re-running with the toolchain present shows the INDEL VCF has **433 and 821
records** in those two genes' CDS±2kb windows, and a precise CDS-exon-overlap screen
(`scripts/check_gene_coding_indels.py`, `results/gene_coding_indel_screen.csv`) finds
**segregating coding indels in the CDS exons of 10/10 target genes** (see below).

Consequences:
- The SNP-only CDS/protein/alignment/variant-table output is **subject to a documented
  caveat for every gene**: coding indels are real and NOT represented in the sequences.
  `screen_indels.sh` was fixed to fail loudly instead of silently reporting 0
  (→ any future gene with a genuine "no indels" claim must pass the loud check).
- This matters for the association step: the **indel genotypes, not just SNP coding
  variants, must be tested against phenotype** — and doing so found the session's
  strongest signal (below).

## Pipeline

1. `scripts/extract_gene_sequences.py` — GFF3 CDS exons + genome FASTA + SNP VCF →
   per-strain spliced CDS FASTA, translated protein FASTA, and a polymorphic-positions
   CSV (invariant CDS columns dropped). No MAFFT: sequences are asserted equal-length
   and built as a plain positional character matrix — safe at the SNP level, with the
   caveat above that real coding indels are NOT captured.
2. `scripts/screen_indels.sh` — read-only screen of the **separate INDEL VCF**
   against each gene's CDS+/-2kb **and now FAILS LOUDLY** if bcftools is unavailable
   or the region returns 0 (fixed 2026-08-26 after the silent-0 bug above).
3. `scripts/check_gene_coding_indels.py` — precise form: counts INDEL-VCF records
   whose POS falls inside a gene's actual CDS exons, and among those how many have a
   non-ref call among the 213 accepted strains. This is the number that distinguishes
   "flank/intergenic noise" from "real coding-indel violation of the SNP-only premise."
4. snpEff (`snpEff/4.3m` module + pre-built `RmucNRRLY2510` database) annotates the SNP
   VCF restricted to each gene's CDS+flank region. **Batch gotchas fixed**: must pass
   `-dataDir <.../snpEff/data>` explicitly (the DB config's `data.dir = ./data/` is
   relative and resolves to the repo root), and use absolute JAVA (module load snpEff
   resets PATH and drops java).
5. `scripts/build_variant_table.py` — parses snpEff `ANN` for the entry matching the
   target gene; joins per-strain genotype, population, phenotype snapshot, lead SNP
   genotype, and Tier A p-value.
6. `scripts/render_alignment_image.py` — static PNG per gene, strains sorted by
   phenotype (pathway genes) or by lead-SNP genotype then phenotype (GWAS-locus genes).
7. **NEW — `scripts/associate_variants_with_phenotype.py`** — population-aware battery
   for every coding-impact (HIGH/MODERATE) SNP variant in each gene × 6 color traits
   (chroma, sat, bright, lab_L, lab_a, lab_b): within-population re-test (with the
   precomputed per-population near-clone collapse maps), fixed-effect inverse-variance
   meta-analysis, population-covariate partial R², BH-FDR across the test matrix,
   and the same `likely_real / likely_population_artifact / ambiguous_underpowered`
   verdict rule as the GWAS disambiguation battery.
8. **NEW — `scripts/associate_indels_with_phenotype.py`** — same battery run on
   segregating in-CDS INDEL genotypes pulled directly from the INDEL VCF (frameshift
   vs in-frame per-allele classification, since these are absent from the SNP tables).

## Pilot genes (carotenoid pathway, corrected) + all 10-gene CDS-indel screen

| gene | product | CDS length | SNP variant sites (region) | polymorphic (panel) | premature stops | CDS-exon indel records (flank±2kb) | segregating in-CDS indels (213 panel) |
|---|---|---|---|---|---|---|---|
| `OM429_003333` | phytoene synthase / lycopene beta-cyclase | 1293 bp | 58 | 49 | 0/213 | 433 | 3 |
| `OM429_003336` | phytoene desaturase (misannotated in GFF3, D-20) | 1929 bp | 81 | 72 | 0/213 | 821 | 84 |
| `OM429_000065` | sat-locus gene (hypothetical) | ... | ... | 157 | — | 585 | 59 |
| `OM429_001521` | lab_b/gwasc (ogg1) | ... | ... | ... | — | 362 | 34 |
| `OM429_001533` | cu_dose/gwas | 1356 bp | 77 | 64 | 0/213 | 499 | 21 |
| `OM429_001415` | lab_L/gwas | ... | ... | ... | — | 397 | 12 |
| `OM429_001430` | lab_L/gwasc | ... | ... | ... | — | 593 | 24 |
| `OM429_003729` | lab_a/gwas (NAP1) | ... | ... | ... | — | 484 | 10 |
| `OM429_002663` | lab_a/gwasc (LYS1) | ... | ... | ... | — | 402 | 2 |
| `OM429_005034` | lab_b/gwas (GYP1) | ... | ... | ... | — | 777 | 92 |

(`...` = populated in the outputs; see `results/*/variant_table.csv`,
`results/gene_coding_indel_screen.csv`.)

**Every one of the 10 target genes has ≥2 segregating indels inside its CDS exons.**
The SNP-only premise does NOT hold for any candidate gene in this panel. Sequences and
variant tables remain the correct description of the *SNP* layer, but any reading of a
gene must check `gene_coding_indel_screen.csv` for the indel layer. Full per-gene
coding-indel site lists are in that CSV's per-gene breakdown (sites printed to stdout).

Reference-strain sanity check (resolves design-doc open question 4): `NRRL_Y-2510` is
homozygous-reference at all 139 CDS SNP variant sites in both pilot genes (0 alt),
confirming the substitution pipeline has no systematic REF/ALT-swap or off-by-one bug.

## Phenotype-association results (adds strict statistical answer to "does any
candidate-gene variation associate with color")

Battery details: within-pop re-test requires ≥3 strains per allele class per population
(after near-clone collapse); meta-analysis requires ≥2 testable populations; FDR across
the full variant×trait matrix. **Reporting convention: "replicated" means
FDR-significant + `likely_real` verdict (directionally consistent across ≥2 testable
populations, meta_p<0.05, partial R²≥0.02).**

### SNP coding variants (954 tests across 10 genes × 6 traits)

**39 FDR-significant, 40 `likely_real`.** Standouts:

| gene | coding variant | trait(s) | meta_p (min) | FDR q | partial R² | replica. |
|---|---|---|---|---|---|---|
| `OM429_000065` (sat locus) | c.839C>T p.Ala280Val (208569, = lead SNP) & c.898C>G p.Leu300Val (208628, same LD) | lab_a, chroma, sat, bright | 8.8e-11 | 1.2e-8 | 0.086–0.10 | 3/3 testable pops, 3/3 same direction |
| `OM429_001521` (ogg1, lab_b/gwasc) | c.735A>C p.Gln245His (503594) | bright | 1.6e-7 | 6.4e-6 | 0.056 | 3/3 |
| `OM429_005034` (GYP1, lab_b/gwas) | c.34T>C p.Tyr12His (608801) | lab_a, chroma, sat | 2.3e-7 | 7.8e-6 | 0.088–0.11 | 3/3 |
| `OM429_002663` (LYS1, lab_a/gwasc) | c.12C>G p.Asn4Lys (884293) | lab_a, chroma, sat | 1.5e-6 | 3.9e-5 | 0.070–0.080 | 3/3 |
| `OM429_001533` (cu_dose/gwas) | c.1280A>T p.Tyr427Phe (544686) | sat, chroma, lab_b | 3.3e-4 | 5.0e-3 | 0.054–0.075 | 3/3 |
| `OM429_000065` | c.2947T>G p.Ser983Ala (211646) | lab_a, bright, chroma | 3.3e-3 | 3.8e-2 | 0.032–0.059 | 3/3 |

**Essentially every significant SNP signal is a coding variant at (or in perfect LD
with) that gene's own GWAS locus-lead region** — they are the *coding layer* of loci
already validated by the population-vs-locus disambiguation (D-17) and MAS gates
(D-18). The new contribution here is pinpointing the specific codon/AA change: e.g.
`OM429_000065` p.Ala280Val and p.Leu300Val are near the GAGCGG-repeat region, and
`OM429_005034` p.Tyr12His sits right at the protein N-terminus.

### Indel genotypes (2046 tests across 10 genes × 6 traits; 1110 frameshift)

**32 FDR-significant. The single strongest candidate-gene association in the whole
analysis is an insertion that the SNP-only pipeline could not even see:**

- **`OM429_000065:208398` — +6 bp (GAGCGG unit) in-frame insertion, meta_p=8.8e-11,
  FDR q=3.6e-8, partial R²=0.086–0.10 for lab_a/chroma/sat, replicated 3/3 testable
  populations, all same direction (insertion → higher chroma/sat/lab_a).** This site
  sits in a coding GAGCGG tandem-repeat, is in **perfect LD (100% concordant) with the
  sat-locus lead SNP scaffold_1:208569** (which is this same gene's p.Ala280Val), and
  is ~170 bp from it inside the same CDS exon. The sat locus therefore has a
  *coding repeat-copy-number polymorphism* as a plausible molecular driver — invisible
  to any SNP-only pipeline. (Signal is repeated-population-consistent; populations 1/2/4
  are near-fixed so only 3 of 6 are testable — report with that caveat.)
- Other FDR-sig indels at `OM429_000065` 208253/208243/208254 (repeat-region cluster),
  `OM429_001521:502667` (+1 bp frameshift, bright), `OM429_003336:175703` (AGTG repeat,
  in-frame ±3 bp; the pilot's only "likely_real" indel, corroborated here).

### Pilot carotenoid genes — what this resolves

- **OM429_003333:169997 splice_donor_variant (HIGH) is MONOMORPHIC in the 213-strain
  panel (0/213 alt)** — it cannot be tested against phenotype and cannot explain
  within-panel color variation (TODO item "splice_donor genotype vs color" is answered:
  no testable variation). The prior flag in `.living/findings/gwas-colocalization.md`
  should be read with this correction in mind.
- **Neither carotenoid gene shows an FDR-significant SNP or indel association with any
  of the 6 color traits** with replication across ≥2 populations, consistent with the
  GWAS port's finding that the carotenoid cluster does NOT co-localize with any
  validated color locus. (Several missense indels at OM429_003336 are directionally
  consistent but ambiguous/underpowered — mostly private-to-population, too few
  testable populations.)

## Open items / honest caveats

- **`screen_indels.sh`'s pilot output was wrong (0 instead of 433/821)**; this session's
  corrected numbers and the new fail-loud behavior supersede the pilot's indel claims
  everywhere. The reach of the SNP-only alignment's "no indels" framing is retracted.
- The near-repeat perfect-LD between the indel and the lead SNP means we cannot separate
  "repeat insertion is causal" from "lead SNP is causal" without denser/functional data —
  but the indel is now a concrete, testable molecular hypothesis (repeat-copy variant in
  a coding region of the sat-locus gene), which is exactly what this pipeline was built to
  surface.
- **13/213 strains lack population assignment** (the 201-strain pop map) and are excluded
  from within-pop tests; the culled-182 panel is the analysis panel, so ~92 carriers in
  the full 213 become ~68 after culling — reported n is per-test.
- `extract_gene_sequences.py`'s 213/213 coverage assertion is present implicitly via the
  strain-list filter (Opus's explicit-assertion recommendation still open).
- **Interactive viewer deferred** (user scoped static PNG+CSV this round).

## Follow-up, not built this round

- Interactive viewer (see CANDIDATE_GENE_ALIGNMENT.md output scoping decision).
- KofamScan results for the 8 GWAS-locus genes (map00906 pathway membership) — KofamScan
  job for the whole genome is pending review; OM429_003333/003336 already confirmed.
- BLASTp/DIAMOND confirmation of medium-confidence pigment-HMM hits before treating any
  as established (`algorithms/functional_annotation` TODO).

## Files

- `scripts/` — extract_gene_sequences.py, screen_indels.sh (fixed), 
  check_gene_coding_indels.py, run_extend_to_gwas_genes.sh (batch driver for the 8
  GWAS-locus genes), build_variant_table.py, render_alignment_image.py,
  associate_variants_with_phenotype.py, associate_indels_with_phenotype.py
- `results/<gene_id>/{dna.fasta, protein.fasta, dna_polymorphic_positions.csv,
  indel_screen.csv, variant_table.csv, variant_table_strain_context.csv, alignment.png}`
  for all 10 genes (8 extended this session: OM429_001533/001415/001430/003729/002663/
  005034/001521/000065)
- `results/gene_coding_indel_screen.csv` — per-gene CDS-exon-exact indel screen
- `results/candidate_gene_phenotype_assoc_all10.csv` (+_within_pop_details.csv),
  `results/candidate_indel_phenotype_assoc_all10.csv` (+_within_pop_details.csv),
  `results/candidate_indel_phenotype_assoc_pilot.csv` — the association step
- `results/PROVENANCE.json` — updated with corrected indel counts, extension genes,
  tool-fix notes
- `results/snpeff_pilot/` — pilot region + snpEff outputs
- `results/extend_batch.log` — the 8-gene extension run log
