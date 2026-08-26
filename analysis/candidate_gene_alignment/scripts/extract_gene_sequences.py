#!/usr/bin/env python3
"""Build per-strain spliced CDS (DNA) and translated protein sequences for a set of
candidate genes, across the 213-strain accepted panel.

Design: docs/superpowers/specs/2026-08-26-candidate-gene-alignment-design.md, reviewed
by an Opus general-review pass and a bioinformatics-expert (fable) consult before
implementation. Key decisions from that review baked in here:

- The genotype VCF is SNP-only (no indels) for every CALLED site, so every strain's CDS
  is the same length as the reference *unless* a real indel exists in this genome that
  was excluded from the SNP VCF entirely -- screened separately by screen_indels.sh, NOT
  handled here (this script will loudly assert equal lengths and fail if that assumption
  is ever violated).
- Missing genotypes are rendered as 'N' in DNA / 'X' in protein (never fall back to the
  reference allele, never silently drop strains) -- both reviewers flagged reference-
  fallback as a real bias risk (masks/erases exactly the variation being characterized).
- Every substituted site is checked against the reference genome base at build time;
  a REF/genome mismatch aborts immediately (loud, not silent) -- catches any
  coordinate-system or scaffold-naming mismatch between the VCF and this GFF3/genome.
- No MAFFT for DNA: sequences are asserted equal-length and built as a plain positional
  matrix (character array), which cannot silently introduce a spurious gap the way a
  general aligner could on data that's already known to be gapless.
"""
import argparse
import gzip
import os
import subprocess
import sys

import pandas as pd

CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def load_cds_exons(gff3_path: str, gene_id: str) -> tuple[str, str, list[tuple[int, int]]]:
    """Returns (scaffold, strand, [(start,end), ...]) in TRANSCRIPT order (i.e. as
    written in the GFF3 -- funannotate writes CDS rows in transcript order: descending
    genomic position for '-' strand genes, ascending for '+')."""
    opener = gzip.open if gff3_path.endswith(".gz") else open
    exons = []
    scaffold, strand = None, None
    with opener(gff3_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "CDS":
                continue
            if f"Parent={gene_id}-T1" not in fields[8]:
                continue
            scaffold = fields[0]
            strand = fields[6]
            exons.append((int(fields[3]), int(fields[4])))
    assert exons, f"no CDS rows found for {gene_id}-T1 in {gff3_path}"
    total_len = sum(e - s + 1 for s, e in exons)
    assert total_len % 3 == 0, f"{gene_id}: total CDS length {total_len} not a multiple of 3"
    return scaffold, strand, exons


def get_genome_seq(samtools_bin: str, genome_fasta: str, scaffold: str, start: int, end: int) -> str:
    r = subprocess.run([samtools_bin, "faidx", genome_fasta, f"{scaffold}:{start}-{end}"],
                        check=True, capture_output=True, text=True)
    lines = r.stdout.strip().splitlines()
    assert lines and lines[0].startswith(">"), f"unexpected samtools faidx output for {scaffold}:{start}-{end}"
    return "".join(lines[1:]).upper()


def get_variants_in_region(bcftools_bin: str, vcf_path: str, scaffold: str, start: int, end: int,
                            strains: list[str]) -> dict[int, dict]:
    """Returns {genomic_pos: {'ref': R, 'alt': A, genotypes: {strain: '0'|'1'|None}}}."""
    fmt = "%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n"
    r = subprocess.run([bcftools_bin, "query", "-r", f"{scaffold}:{start}-{end}", "-f", fmt, vcf_path],
                        check=True, capture_output=True, text=True)
    variants = {}
    strain_set = set(strains)
    for line in r.stdout.splitlines():
        parts = line.split("\t")
        pos, ref, alt = int(parts[1]), parts[2], parts[3]
        if "," in alt:
            continue  # multiallelic -- skip for this simple biallelic substitution model
        genos = {}
        for tok in parts[4:]:
            if "=" not in tok:
                continue
            s, gt = tok.split("=")
            if s not in strain_set:
                continue
            gt = gt.replace("|", "/").split("/")[0]
            genos[s] = gt if gt in ("0", "1") else None
        variants[pos] = dict(ref=ref, alt=alt, genotypes=genos)
    return variants


def build_gene_sequences(gene_id: str, scaffold: str, strand: str, exons: list[tuple[int, int]],
                          genome_fasta: str, vcf_path: str, strains: list[str],
                          samtools_bin: str, bcftools_bin: str) -> dict:
    exon_lo = min(s for s, e in exons)
    exon_hi = max(e for s, e in exons)
    variants = get_variants_in_region(bcftools_bin, vcf_path, scaffold, exon_lo, exon_hi, strains)

    # Load reference sequence per exon (genomic forward orientation) and validate REF alleles.
    exon_ref_seqs = []
    for (s, e) in exons:
        ref_seq = get_genome_seq(samtools_bin, genome_fasta, scaffold, s, e)
        exon_ref_seqs.append(list(ref_seq))
        for pos, v in variants.items():
            if s <= pos <= e:
                genome_base = ref_seq[pos - s]
                assert genome_base == v["ref"], (
                    f"{gene_id}: VCF REF={v['ref']} at {scaffold}:{pos} does not match "
                    f"genome base {genome_base} -- coordinate system or scaffold mismatch, aborting"
                )
    n_variant_sites = sum(1 for pos in variants if exon_lo <= pos <= exon_hi)

    per_strain_cds = {}
    for strain in strains:
        exon_strain_seqs = []
        for (s, e), ref_list in zip(exons, exon_ref_seqs):
            strain_exon = list(ref_list)  # copy
            for pos, v in variants.items():
                if not (s <= pos <= e):
                    continue
                gt = v["genotypes"].get(strain)
                offset = pos - s
                if gt is None:
                    strain_exon[offset] = "N"
                elif gt == "0":
                    pass  # already reference
                elif gt == "1":
                    strain_exon[offset] = v["alt"]
            exon_seq = "".join(strain_exon)
            if strand == "-":
                exon_seq = revcomp(exon_seq)
            exon_strain_seqs.append(exon_seq)
        per_strain_cds[strain] = "".join(exon_strain_seqs)

    lengths = set(len(v) for v in per_strain_cds.values())
    assert len(lengths) == 1, f"{gene_id}: strains have unequal CDS lengths {lengths} -- unexpected given SNP-only input"

    return dict(gene_id=gene_id, scaffold=scaffold, strand=strand, exons=exons,
                cds_length=lengths.pop(), n_variant_sites=n_variant_sites,
                variants=variants, sequences=per_strain_cds)


def translate(cds: str) -> tuple[str, bool, int | None]:
    """Returns (protein, premature_stop, stop_codon_index). Codons containing 'N' -> 'X'."""
    protein = []
    n_codons = len(cds) // 3
    premature_stop = False
    stop_idx = None
    for i in range(n_codons):
        codon = cds[i * 3:(i + 1) * 3]
        if "N" in codon:
            protein.append("X")
            continue
        aa = CODON_TABLE.get(codon, "X")
        if aa == "*":
            stop_idx = i
            if i != n_codons - 1:
                premature_stop = True
            break
        protein.append(aa)
    return "".join(protein), premature_stop, stop_idx


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene-ids", nargs="+", required=True)
    ap.add_argument("--gff3", required=True)
    ap.add_argument("--genome-fasta", required=True)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--strain-list", required=True, help="accepted VCF sample IDs, one per line")
    ap.add_argument("--samtools", default="samtools")
    ap.add_argument("--bcftools", default="bcftools")
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    strains = sorted(set(l.strip() for l in open(args.strain_list) if l.strip()))
    print(f"Strains: {len(strains)}")

    for gene_id in args.gene_ids:
        print(f"\n=== {gene_id} ===")
        scaffold, strand, exons = load_cds_exons(args.gff3, gene_id)
        print(f"  {scaffold} {strand} {len(exons)} exon(s)")
        result = build_gene_sequences(gene_id, scaffold, strand, exons, args.genome_fasta,
                                       args.vcf, strains, args.samtools, args.bcftools)
        print(f"  CDS length: {result['cds_length']} bp, {result['n_variant_sites']} variant sites")

        gene_dir = os.path.join(args.out_dir, gene_id)
        os.makedirs(gene_dir, exist_ok=True)

        with open(os.path.join(gene_dir, "dna.fasta"), "w") as f:
            for s in strains:
                f.write(f">{s}\n{result['sequences'][s]}\n")

        n_missing_by_strain = {s: result["sequences"][s].count("N") for s in strains}
        n_premature_stop = 0
        with open(os.path.join(gene_dir, "protein.fasta"), "w") as f:
            protein_lengths = set()
            for s in strains:
                protein, premature, stop_idx = translate(result["sequences"][s])
                protein_lengths.add(len(protein))
                if premature:
                    n_premature_stop += 1
                    f.write(f">{s} PREMATURE_STOP_at_codon_{stop_idx}\n{protein}\n")
                else:
                    f.write(f">{s}\n{protein}\n")
        print(f"  premature stop codons: {n_premature_stop}/{len(strains)} strains")
        print(f"  protein lengths seen: {sorted(protein_lengths)}")
        print(f"  strains with >=1 missing (N) site: {sum(1 for v in n_missing_by_strain.values() if v > 0)}")

        # Positional DNA matrix (strain x CDS position), only for polymorphic positions
        # (dropping invariant columns keeps this readable -- 1500+ constant columns add
        # nothing).
        poly_positions = [i for i in range(result["cds_length"])
                          if len(set(result["sequences"][s][i] for s in strains)) > 1]
        matrix_rows = []
        for s in strains:
            row = {"strain": s}
            for i in poly_positions:
                row[f"cds_pos_{i+1}"] = result["sequences"][s][i]
            matrix_rows.append(row)
        pd.DataFrame(matrix_rows).to_csv(os.path.join(gene_dir, "dna_polymorphic_positions.csv"), index=False)
        print(f"  polymorphic CDS positions (0 constant columns dropped): {len(poly_positions)}")


if __name__ == "__main__":
    main()
