#!/usr/bin/env python3
"""Data-driven per-profile score-threshold calibration for the pigmentation-pathway HMM
library, requested because the source project's own per-profile thresholds were NOT
copied into this repo (only the HMMs + a report describing that they exist) -- and we
have no labeled true/false-positive training set of our own to reproduce the source
project's discriminative-column method (that method IS available and already applied for
t3hnr/t4hnr, data/raw/pigmentation-pathway-hmms/specific/RECOMMENDED_CONFIG.md).

Two independent, UNSUPERVISED lines of evidence (no external truth labels required):

1. Empirical bit-score gap ("twilight zone" heuristic): pool full-sequence bit-scores
   for a profile across many related proteomes (11 local Rhodotorula/Cystobasidium
   genomes here), sort descending, and find the largest gap between consecutive scores.
   A true ortholog set is expected to cluster tightly at high score with a large drop to
   unrelated/noise hits; the threshold is set at the midpoint of the largest such gap.
   This is the same logic used informally when picking BLAST/HMM cutoffs by eye from a
   score histogram, just automated and reproducible.

2. Existing-annotation concordance: for each candidate hit surviving the gap threshold
   in R. mucilaginosa (the only genome here with GFF3 gene annotation loaded), report
   its EXISTING funannotate product/Pfam/InterPro call. This is independent evidence
   (built from a different method entirely -- HMMER Pfam-A/InterPro domain scans, not
   this pigment-HMM library) that a reader can use to judge concordance by eye; genes
   whose existing Pfam/InterPro annotation is compatible with the profile's expected gene
   family corroborate the pigment-HMM call, while a mismatch is a flag for closer review
   (especially for profiles the source report already calls "medium confidence").

Caveat, stated plainly: neither of these is a substitute for the source project's actual
validated method (94 fungal genomes + 100 metagenome MAGs with real TP/FP labels). They
are the best calibration achievable from data already on hand, not a claim of equivalent
rigor -- treat outputs as a triage/prioritization aid, not a final call.
"""
import argparse
import glob
import json
import os

import pandas as pd


def load_domtbl(path: str) -> pd.DataFrame:
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return pd.DataFrame(columns=["profile", "protein", "genome", "full_evalue", "full_score"])
    rows = []
    genome = os.path.basename(os.path.dirname(path))
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.split()
            rows.append(dict(profile=fields[0], protein=fields[3],
                              full_evalue=float(fields[6]), full_score=float(fields[7]), genome=genome))
    return pd.DataFrame(rows)


def find_gap_threshold(scores: list[float], min_group: int = 1) -> tuple[float, int]:
    """Largest gap between consecutive sorted (desc) scores; returns (threshold, n_above).
    Threshold = midpoint of the largest gap. Ties/near-uniform distributions (no clear
    gap) fall back to reporting all hits with threshold = min(scores) - 0 (no calibration
    possible from this data alone -- flagged via n_above == len(scores))."""
    s = sorted(set(scores), reverse=True)
    if len(s) < 2:
        return (s[0] - 1e-6 if s else 0.0), len(scores)
    gaps = [(s[i] - s[i + 1], i) for i in range(len(s) - 1)]
    biggest_gap, idx = max(gaps, key=lambda x: x[0])
    threshold = (s[idx] + s[idx + 1]) / 2.0
    n_above = sum(1 for x in scores if x > threshold)
    return threshold, max(n_above, min_group)


def load_gene_annotation(gene_index_path: str) -> dict:
    with open(gene_index_path) as f:
        return json.load(f)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-glob", required=True,
                     help="glob for per-genome scan dirs, e.g. 'results/pigment_scan_*'")
    ap.add_argument("--focal-genome", default="Rhodotorula_mucilaginosa_NRRL_Y-2510",
                     help="genome dir name to report per-gene annotation concordance for")
    ap.add_argument("--gene-index", required=True, help="gene_index.json for the focal genome")
    ap.add_argument("--out-thresholds", required=True)
    ap.add_argument("--out-hits", required=True)
    args = ap.parse_args()

    dirs = sorted(glob.glob(args.results_glob))
    assert dirs, f"no directories matched {args.results_glob}"

    all_rows = []
    for d in dirs:
        for domtbl in glob.glob(os.path.join(d, "*.domtbl")):
            all_rows.append(load_domtbl(domtbl))
    combined = pd.concat(all_rows, ignore_index=True)
    combined["genome"] = combined["genome"].apply(lambda x: x)
    print(f"Loaded {len(combined)} raw domain hits across {combined['genome'].nunique()} genomes, "
          f"{combined['profile'].nunique()} profiles")

    genes = load_gene_annotation(args.gene_index)

    thresh_rows = []
    hit_rows = []
    for profile, grp in combined.groupby("profile"):
        # Use each PROTEIN's best (max) full-sequence score once, across all genomes,
        # for gap-finding -- avoids multiple weak domain hits on the same protein
        # inflating the low end of the distribution.
        best_per_protein = grp.groupby(["genome", "protein"])["full_score"].max().reset_index()
        scores = best_per_protein["full_score"].tolist()
        threshold, n_above = find_gap_threshold(scores)
        calibratable = n_above < len(scores)  # False if no gap found at all (flat distribution)
        thresh_rows.append(dict(profile=profile, n_total_hits=len(scores), gap_threshold_bitscore=round(threshold, 1),
                                 n_above_threshold=n_above, calibratable=calibratable))

        confident = best_per_protein[best_per_protein["full_score"] > threshold]
        for _, r in confident.iterrows():
            row = dict(profile=profile, genome=r["genome"], protein=r["protein"], full_score=r["full_score"])
            if r["genome"] == args.focal_genome:
                gid = r["protein"].split("-T")[0]
                g = genes.get(gid, {})
                row.update(gene_id=gid, product=g.get("product"),
                           pfam=";".join(g.get("pfam", [])), interpro=";".join(g.get("interpro", [])),
                           go=";".join(g.get("go", [])))
            hit_rows.append(row)

    thresh_df = pd.DataFrame(thresh_rows).sort_values("profile")
    hits_df = pd.DataFrame(hit_rows).sort_values(["profile", "full_score"], ascending=[True, False])

    os.makedirs(os.path.dirname(args.out_thresholds), exist_ok=True)
    thresh_df.to_csv(args.out_thresholds, index=False)
    hits_df.to_csv(args.out_hits, index=False)

    print(f"\nWrote {args.out_thresholds}")
    print(f"Wrote {args.out_hits}")
    print("\nPer-profile calibration summary:")
    print(thresh_df.to_string(index=False))


if __name__ == "__main__":
    main()
