# Algorithms Manifest

<!-- Add entries below using the appropriate manifest entry template. -->

### functional-annotation
```yaml
name: functional-annotation
question: How can genes in a funannotate-style genome (starting with R. mucilaginosa NRRL Y-2510) be functionally annotated for a specific pathway of interest (pigmentation) and, as a reusable baseline, genome-wide KEGG Orthology -- when a pre-built HMM library exists but its calibrated per-profile thresholds were not copied over?
inputs:
  - proteome FASTA (funannotate-style, e.g. Rhodotorula_mucilaginosa_NRRL_Y-2510.proteins.fa)
  - data/raw/pigmentation-pathway-hmms/ (28-profile HMM library, imported resource)
  - gene_index.json (product/GO/InterPro/Pfam per gene)
  - /bigdata/stajichlab/shared/lib/gtotree/kofamscan_data/ (KOfam HMM DB, already on shared storage)
scripts:
  - algorithms/functional_annotation/scripts/run_pigment_hmm_scan.sh
  - algorithms/functional_annotation/scripts/calibrate_pigment_hmm_thresholds.py
  - algorithms/functional_annotation/scripts/run_kofamscan.sh
outputs:
  - algorithms/functional_annotation/results/pigment_scan_<genome>/*.domtbl
  - algorithms/functional_annotation/results/calibrated_thresholds.csv
  - algorithms/functional_annotation/results/confident_hits.csv
  - algorithms/functional_annotation/results/kofamscan/<genome>.{detail,mapper}.txt
reproduce: see algorithms/functional_annotation/FUNCTIONAL_ANNOTATION.md "Usage"
status: in-progress (pigment-pathway HMM scan + calibration complete across 11 local Rhodotorula/Cystobasidium proteomes; genome-wide KofamScan for R. mucilaginosa submitted, pending completion)
key_findings:
  - Product-name keyword search (the first approach tried, analysis/gwas/scripts/finemap_candidate_genes.py) found NO carotenoid-pathway gene near any GWAS locus -- but direct HMM homology search immediately found the actual fungal carotenoid gene cluster on scaffold_7:168525-175999, which keyword search had missed entirely because both genes are misannotated in the existing GFF3 (OM429_003333 = "hypothetical protein" but scores 689.9/751.1 against crt_fungal_psy/crt_fungal_lcy -- the classic bifunctional phytoene synthase/lycopene cyclase fusion; OM429_003336 = "transcription factor, PHD finger motif" but scores 877.4 against crt_fungal_pds, strongly suggesting it is actually the phytoene desaturase).
  - This carotenoid gene cluster is NOT colocalized with any of the 8 GWAS-validated likely_real loci (D-17/D-18); the closest genome-wide signal is lab_L's scaffold_7:172154 (p=1.1e-7, FDR-significant, rank 104/496,358) sitting in the intergenic region between the two genes -- but this SNP FAILS the within-population validation test (only 1 of 6 populations had testable allelic variance; population 4 is ~98% fixed for the allele), so it is reported as a tantalizing but statistically unconfirmed lead, not a validated finding.
  - Data-driven bit-score-gap calibration (no labeled training data available) found large, clean score gaps (600+ bits) for the three fungal carotenoid profiles across 11 genomes -- high-confidence single-copy orthologs each -- while several other profiles (mysB, single hits only) could not be meaningfully calibrated from data on hand.
tags: [functional-annotation, hmmer, hmmscan, kofamscan, kegg, pigmentation, carotenoid, gene-annotation, calibration, reusable, algorithm]
```
