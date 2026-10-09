# RNA-seq SRA scoping — user-supplied BioProjects

**Status**: scoping only (no expression analysis run). User supplied 3 candidate BioProjects
directly (2026-08-27) rather than waiting on the parallel background research-agent scan
(dispatched via the `research` skill for a broader NCBI SRA/GEO search; if/when it completes
separately, reconcile its list against this one rather than replacing it).

**Purpose**: identify public R. mucilaginosa RNA-seq datasets with ≥3 biological replicates
per condition, as future orthogonal evidence for GWAS/fine-mapping candidate genes
currently annotated only as "hypothetical protein" (e.g. `OM429_005716`/`OM429_005715`
near the `scaffold_13:810026` locus, see `GWAS.md` §20/§23; also the color/copper
candidates in `analysis/candidate_gene_alignment/`).

**Method**: pulled BioProject, SRA experiment, and BioSample metadata directly from NCBI
(bioproject/sra/biosample web records) via WebFetch for each accession the user supplied.
No sequence data downloaded, no alignment or expression analysis performed.

## Summary table

| BioProject | Strain | Design | Reps/group | Meets ≥3? | Reference-strain match |
|---|---|---|---|---|---|
| PRJNA628936 | unspecified (biocontrol isolate) | control vs. chitosan-treated | **2** | No | No (NRRL Y-2510 not used) |
| PRJNA954140 | BMU419 | 0h vs. 24h cultivation timepoint | **3** | Yes | No |
| PRJNA1451956 | LWJJ06 (citrinin-degrading, fermented-tea isolate) | control vs. 10 vs. 50 [conc. units] citrinin | **4** | Yes, best-replicated | No |

**None of the three use this project's reference strain (NRRL Y-2510).** Any expression
comparison against the `OM429_*` gene models used throughout this project's GWAS work
would require an orthology-mapping step (e.g. reciprocal-best-BLAST or OrthoFinder against
each study's own assembly/annotation, if one exists) before coordinate-based expression
counts could be trusted — this was not checked (would require knowing whether BMU419 and
LWJJ06 have public genome assemblies, which was out of scope for this pass).

## PRJNA628936 — chitosan-induced biocontrol transcriptome

- **Title**: "Transcriptome analysis reveals the mechanisms involved in the enhanced
  biocontrol efficacy of Rhodotorula mucilaginosa induced by chitosan"
- **Submitter**: Jiangsu University; registered 2020-04-28
- **Design**: 4 SRA experiments = 2 conditions × 2 biological replicates (library names
  `YA-1`/`YA-2` = control in nutrient yeast dextrose broth, `YB-1`/`YB-2` = chitosan-treated)
- **Samples**: SRX8189311 (YA-1, SRR11625205), SRX8189312 (YA-2, SRR11625204),
  SRX8189313 (YB-1, SRR11625203), SRX8189314 (YB-2, SRR11625202)
- **Platform**: Illumina HiSeq X Ten, paired-end, ~28-31M spots/sample
- **No linked publication found** on the BioProject/SRA pages themselves.
- **Verdict**: **Fails the ≥3-replicate bar (n=2/group).** Not useful for a confident
  differential-expression test; lowest priority of the three, and least thematically
  relevant (biocontrol/chitosan induction, not stress/toxin dose-response or growth
  timepoint).

## PRJNA954140 — R. mucilaginosa BMU419 growth-timepoint transcriptome

- **Title**: "Rhodotorula mucilaginosa strain:BMU419 Raw sequence reads"
- **Submitter**: Binzhou Medical University; registered 2023-04-10
- **Design**: 6 SRA experiments = 2 timepoints × 3 biological replicates. Library names
  P01/P02/P03 = "0 h of cultivation time"; P241/P242/P243 = "24 h of cultivation time".
- **Samples**: SRX19924329 (P01, SAMN34134732), SRX19924330 (P02, SAMN34134733),
  SRX19924332 (P241, SAMN34134735) confirmed directly; P03/P242/P243 (SRX19924331,
  SRX19924333, SRX19924334) inferred from the same naming series, not individually
  re-verified.
- **Platform**: Illumina HiSeq 2000, paired-end, polyA selection.
- **No linked publication identified** on the pages fetched.
- **Verdict**: **Meets the ≥3-replicate bar.** Most thematically relevant of the three to
  this project's growth-rate/resilience traits (`resilience_30`/`AUC_30`,
  `scaffold_13:810026`) since it's a direct early-vs-late growth-phase comparison — a
  natural place to check whether an ortholog of `OM429_005716` or `OM429_005715` (`ubc12`)
  changes expression over the growth curve, if an ortholog can be identified.

## PRJNA1451956 — R. mucilaginosa LWJJ06 citrinin dose-response transcriptome

- **Title**: "Transcriptome analysis of citrinin-degrading Rhodotorula mucilaginosa LWJJ06"
- **Submitter**: Guilin University of Technology; registered 2026-04-10. Strain deposited
  at China General Microbiological Culture Collection Center as CGMCC No. 35099,
  isolated from fermented tea.
- **Design**: 12 SRA experiments = 3 conditions × 4 biological replicates, confirmed by
  checking every one of the 12 SRX library names directly:
  - Control: `CK1` (SRX32904507), `CK2` (SRX32904508), `CK3` (SRX32904511), `CK4` (SRX32904512)
  - 10 [conc. units] citrinin: `CIT10_1` (SRX32904513), `CIT10_2` (SRX32904514),
    `CIT10_3` (SRX32904515), `CIT10_4` (SRX32904516)
  - 50 [conc. units] citrinin: `CIT50_1` (SRX32904517), `CIT50_2` (SRX32904518),
    `CIT50_3` (SRX32904509), `CIT50_4` (SRX32904510)
  - (Exact concentration units, e.g. mg/L vs µg/mL, not stated on the pages fetched —
    check the BioSample records or an associated publication before citing a number.)
- **Platform**: Illumina HiSeq 2500, paired-end cDNA libraries.
- **No linked publication identified** on the pages fetched (very recently registered,
  2026-04; a paper may not yet be public).
- **Verdict**: **Best-replicated of the three (n=4/group, exceeds the ≥3 bar).** A clean
  two-dose xenobiotic-stress dose-response design (control/low/high), methodologically
  analogous in spirit to this project's own copper dose-response phenotyping even though
  the specific compound (a mycotoxin, not a heavy metal) differs. Worth checking whether
  citrinin and copper stress-response pathways overlap (e.g. general oxidative-stress/
  detoxification genes) before assuming direct relevance to the copper-specific
  candidates; more directly useful as a general stress-responsiveness check for any
  candidate gene (does it respond to xenobiotic stress at all?) than as a copper-specific
  validation.

## Recommendation

If pursuing expression-based follow-up on the `scaffold_13:810026` candidates
(`OM429_005716`, `OM429_005715`/`ubc12`) or the color/copper candidates: start with
**PRJNA1451956** (best replication, n=4/group, clear dose-response structure) for a
general stress-responsiveness check, and **PRJNA954140** (n=3/group, growth timepoint)
for a growth-phase-expression check specifically relevant to the resilience/AUC traits.
**PRJNA628936** is under-replicated (n=2) and should be treated as supplementary evidence
at best, not a standalone test. All three require an orthology-mapping step first since
none uses the NRRL Y-2510 reference strain — this has not been done and is the next
concrete task if this line of evidence is pursued.
