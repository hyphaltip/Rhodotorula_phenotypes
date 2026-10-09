#!/usr/bin/env python3
"""Per-strain heterozygosity + depth cross-check for haploid ploidy validation.

Usage:
  pixi run python3 analysis/gwas/scripts/check_ploidy.py \
      --strains analysis/gwas/results/strain_reconciliation/strain_match_table.reviewed.csv \
      --vcf data/raw/genotypes/RmucY2510_v2/RmucY2510_v2.All.SNP.combined_selected.vcf.gz \
      --mosdepth-dir /bigdata/stajichlab/shared/projects/Population_Genomics/Rhodotorula_mucilaginosa_NRRLY2510/coverage/mosdepth \
      --cram-dir /bigdata/stajichlab/shared/projects/Population_Genomics/Rhodotorula_mucilaginosa_NRRLY2510/aln \
      --out analysis/gwas/results/ploidy_check/ploidy_flags.csv

Thresholds are panel-relative, not fixed (spec S3, quant-genetics review fix #2):
  - het outlier: het_rate > median(het_rate) + 3 * MAD(het_rate)
  - diploidy depth band: depth_ratio in [1.5, 2.5]
  - contamination band: depth_ratio > 2.5
Minimum callable sites: 100,000 (below this, verdict is 'insufficient_data' rather
than a het-rate call — depth-starved samples give noisy het estimates).
"""
import argparse
import csv
import os
import statistics
import subprocess

MIN_CALLABLE_SITES = 100_000
DIPLOID_BAND = (1.5, 2.5)


def accepted_strains(path: str) -> list[str]:
    ids = []
    with open(path) as f:
        for row in csv.DictReader(f):
            tier = row.get("tier", "")
            decision = row.get("decision", "")
            if tier in ("exact", "normalized") or decision == "accept":
                if row["vcf_sample_id"]:
                    ids.append(row["vcf_sample_id"])
    assert len(ids) > 0, f"no accepted strains found in {path}"
    return sorted(set(ids))


def het_and_callable_all(vcf_path: str, samples: list[str]) -> dict[str, tuple[int, int]]:
    """Het count and non-missing (callable) genotype count for every requested sample,
    via a single bcftools stats -s - pass over all 422 VCF samples (one bcftools stats
    call per sample took ~51s; a single -s - pass over all samples takes ~56s total).

    This VCF encodes genotypes haploid (GT like "0"/"1", not diploid "0/1") for a
    haploid organism, so per-sample callable sites live in nHapRef+nHapAlt, not
    nRefHom+nNonRefHom (those stay 0 for haploid-encoded samples); a genuinely
    diploid/aneuploid sample would show up via nHets > 0 instead.
    """
    out = subprocess.run(
        ["bcftools", "stats", "-s", "-", vcf_path],
        check=True, capture_output=True, text=True,
    )
    wanted = set(samples)
    result = {}
    for line in out.stdout.splitlines():
        if not line.startswith("PSC\t"):
            continue
        fields = line.split("\t")
        # PSC: [0]PSC [1]id [2]sample [3]nRefHom [4]nNonRefHom [5]nHets [6]nTransitions
        #      [7]nTransversions [8]nIndels [9]avgDepth [10]nSingletons [11]nHapRef
        #      [12]nHapAlt [13]nMissing
        sample = fields[2]
        if sample not in wanted:
            continue
        n_ref_hom = int(fields[3])
        n_nonref_hom = int(fields[4])
        n_hets = int(fields[5])
        n_hap_ref = int(fields[11])
        n_hap_alt = int(fields[12])
        het = n_hets
        callable_sites = n_ref_hom + n_nonref_hom + n_hets + n_hap_ref + n_hap_alt
        result[sample] = (het, callable_sites)
    missing_samples = wanted - result.keys()
    assert not missing_samples, f"bcftools stats returned no PSC row for samples: {sorted(missing_samples)}"
    return result


def mosdepth_mean(mosdepth_dir: str, cram_dir: str, sample: str) -> float:
    summary = os.path.join(mosdepth_dir, f"{sample}.10000bp.mosdepth.summary.txt")
    if not os.path.exists(summary):
        raise FileNotFoundError(
            f"mosdepth summary missing for {sample}: {summary}. "
            f"Generate it from {os.path.join(cram_dir, sample + '.cram')} with "
            f"`mosdepth --by 10000 {sample} <cram>` before rerunning this script."
        )
    with open(summary) as f:
        for line in f:
            fields = line.split("\t")
            if fields[0] == "total":
                return float(fields[3])  # mean depth column
    raise ValueError(f"no 'total' row in {summary}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--strains", required=True)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--mosdepth-dir", required=True)
    ap.add_argument("--cram-dir", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    strains = accepted_strains(args.strains)

    het_callable = het_and_callable_all(args.vcf, strains)

    records = []
    for s in strains:
        het, callable_sites = het_callable[s]
        try:
            depth = mosdepth_mean(args.mosdepth_dir, args.cram_dir, s)
        except FileNotFoundError as e:
            print(f"WARNING: {e}")
            depth = None
        het_rate = het / callable_sites if callable_sites > 0 else None
        records.append({"strain_id": s, "het_rate": het_rate, "callable_sites": callable_sites, "mean_depth": depth})

    valid_het = [r["het_rate"] for r in records if r["het_rate"] is not None and r["callable_sites"] >= MIN_CALLABLE_SITES]
    valid_depth = [r["mean_depth"] for r in records if r["mean_depth"] is not None]
    assert len(valid_het) >= 10, "fewer than 10 strains have usable het data — panel-relative thresholds are unreliable"
    assert len(valid_depth) >= 10, "fewer than 10 strains have depth data — panel-relative thresholds are unreliable"

    median_het = statistics.median(valid_het)
    mad_het = statistics.median([abs(x - median_het) for x in valid_het])
    median_depth = statistics.median(valid_depth)
    het_outlier_cutoff = median_het + 3 * mad_het

    out_rows = []
    for r in records:
        if r["callable_sites"] < MIN_CALLABLE_SITES or r["mean_depth"] is None:
            flag = "insufficient_data"
            depth_ratio = None
        else:
            depth_ratio = r["mean_depth"] / median_depth
            het_outlier = r["het_rate"] > het_outlier_cutoff
            in_diploid_band = DIPLOID_BAND[0] <= depth_ratio <= DIPLOID_BAND[1]
            if depth_ratio > DIPLOID_BAND[1]:
                flag = "watch_contamination"
            elif het_outlier and in_diploid_band:
                flag = "likely_diploid_or_mixed"
            elif het_outlier or in_diploid_band:
                flag = "watch"
            else:
                flag = "haploid_ok"
        out_rows.append({
            "strain_id": r["strain_id"],
            "het_rate": r["het_rate"],
            "callable_sites": r["callable_sites"],
            "mean_depth": r["mean_depth"],
            "depth_ratio_to_panel_median": depth_ratio,
            "flag": flag,
            "panel_median_het": median_het,
            "panel_mad_het": mad_het,
            "panel_median_depth": median_depth,
        })

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()))
        writer.writeheader()
        writer.writerows(out_rows)

    from collections import Counter
    counts = Counter(r["flag"] for r in out_rows)
    print(f"Strains checked: {len(out_rows)}")
    for flag, n in counts.items():
        print(f"  {flag}: {n}")
    print(f"Panel median het: {median_het:.6f}  MAD: {mad_het:.6f}  cutoff: {het_outlier_cutoff:.6f}")
    print(f"Panel median depth: {median_depth:.2f}")
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
