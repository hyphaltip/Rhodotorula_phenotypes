# Discriminative HMM Configuration Summary

## Problem
t4hnr and t3hnr belong to the SDR superfamily (PF00106), one of the largest
enzyme families across all domains of life. Their original HMMs match general
bacterial SDR enzymes, producing ~55% false positive rate in desert metagenome
MAGs.

## Solution: Discriminative Column Selection

Identified alignment columns where DHN-melanin reductases (fungal seed
alignment) have conserved residues that differ from bacterial SDR false
positives. Built new HMMs from only these discriminating columns plus
flanking context.

### Method
1. Extracted 200 false-positive SDR sequences from desert MAGs
2. Aligned them to the seed HMM coordinate frame using hmmalign --mapali
3. Filtered to match-only columns (331 for t4hnr, 342 for t3hnr)
4. Identified discriminating columns (conserved in seeds, different in FPs)
5. Built HMMs from top N discriminating columns + flanking residues
6. Optimized E-value threshold for sensitivity/specificity balance

## Recommended Configurations

### t4hnr (tetrahydroxynaphthalene reductase)
- HMM: t4hnr_disc_top15_f3.hmm (15 discriminating columns + 3 flanking = 74 cols)
- E-value: 1.0
- Results (94 fungal genomes, 100 desert MAGs):
  - True positive rate: 80.9% (76/94)
  - False positive rate: 2.0% (2/100)
  - FP ratio: 0.02 (vs 0.56 original — 28x improvement)
- Discriminating columns: 34 found at conservation >= 0.80

### t3hnr (trihydroxynaphthalene reductase)
- HMM: t3hnr_disc_top20_f4.hmm (20 discriminating columns + 4 flanking = ~100 cols)
- E-value: 1e-5
- Results (94 fungal genomes, 100 desert MAGs):
  - True positive rate: 85.1% (80/94)
  - False positive rate: 10.0% (10/100)
  - FP ratio: 0.12 (vs 0.57 original — 5x improvement)
- Discriminating columns: 74 found at conservation >= 0.80

## Comparison with Other DHN-melanin Genes

| Gene | Original FP ratio | Discriminative FP ratio | Improvement |
|------|-------------------|-------------------------|-------------|
| scd  | 0.00 (gold std)   | N/A (already specific)  | —           |
| ayg1 | 0.03              | N/A (already specific)  | —           |
| pks_melanin | 0.05       | N/A (already specific)  | —           |
| t4hnr | 0.56             | 0.02                    | 28x         |
| t3hnr | 0.57             | 0.12                    | 5x          |

## Limitations
- t4hnr: 19% of fungal true positives not detected (divergent substrate-binding)
- t3hnr: 15% of fungal true positives not detected
- For maximum sensitivity, combine with pathway context filter:
  Count t4hnr/t3hnr as DHN-melanin if EITHER:
  (a) Detected by discriminative HMM, OR
  (b) Detected by original HMM AND scd or pks_melanin also present in same genome/MAG

## Files
- HMMs: /mnt/results/hmms/specific/t4hnr_disc_top15_f3.hmm, t3hnr_disc_top20_f4.hmm
- Alignments: /mnt/results/hmms/specific/t4hnr_disc_top15_f3.fasta, t3hnr_disc_top20_f4.fasta
- Discriminating columns: /mnt/results/hmms/specific/t4hnr_discriminating_columns_v2.tsv, t3hnr_discriminating_columns_v2.tsv
- Full results: /mnt/results/hmms/specific/full_validation_results.tsv
