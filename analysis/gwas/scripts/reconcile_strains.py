#!/usr/bin/env python3
"""Reconcile phenotype-metadata strain IDs against genotype VCF sample IDs.

Usage:
  pixi run python3 analysis/gwas/scripts/reconcile_strains.py \
      --pheno data/metadata/Copper.Strain_info.csv \
      --pheno-col <STRAIN_ID_COLUMN_FROM_STEP_1> \
      --vcf data/raw/genotypes/RmucY2510_v2/RmucY2510_v2.All.SNP.combined_selected.vcf.gz \
      --out-dir analysis/gwas/results/strain_reconciliation

Never auto-accepts a fuzzy match or auto-drops an unmatched strain (spec S2):
fuzzy/unmatched rows are written to NEEDS_REVIEW.md for a human to adjudicate in
strain_match_table.reviewed.csv before any GWAS step uses this strain list.
"""
import argparse
import csv
import re
import subprocess
import sys
from collections import defaultdict
from difflib import SequenceMatcher


def normalize(s: str) -> str:
    s = s.strip().lower()
    s = re.sub(r"[_\-\s]+", "_", s)
    s = re.sub(r"^0+(?=\d)", "", s)  # strip leading zeros before a digit run
    return s


def vcf_sample_ids(vcf_path: str) -> list[str]:
    out = subprocess.run(
        ["bcftools", "query", "-l", vcf_path],
        check=True, capture_output=True, text=True,
    )
    ids = [line.strip() for line in out.stdout.splitlines() if line.strip()]
    assert len(ids) > 0, f"bcftools returned 0 sample IDs from {vcf_path}"
    return ids


def suffix_match(pheno_id: str, vcf_ids: list[str]) -> str | None:
    """Longest vcf_id that is a suffix of pheno_id or vice versa (generalizes the
    idea_09_phylogeny.R tip<->strain_code suffix trick)."""
    candidates = []
    for vid in vcf_ids:
        if pheno_id.endswith(vid) or vid.endswith(pheno_id):
            candidates.append(vid)
    if not candidates:
        return None
    return max(candidates, key=len)


def fuzzy_match(pheno_id: str, vcf_ids: list[str], threshold: float = 0.85):
    best_id, best_score = None, 0.0
    for vid in vcf_ids:
        score = SequenceMatcher(None, normalize(pheno_id), normalize(vid)).ratio()
        if score > best_score:
            best_id, best_score = vid, score
    if best_score >= threshold:
        return best_id, best_score
    return None, best_score


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pheno", required=True)
    ap.add_argument("--pheno-col", required=True)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    with open(args.pheno) as f:
        reader = csv.DictReader(f)
        assert args.pheno_col in reader.fieldnames, (
            f"--pheno-col {args.pheno_col!r} not found in {args.pheno} header: "
            f"{reader.fieldnames}"
        )
        pheno_ids = sorted({row[args.pheno_col].strip() for row in reader if row[args.pheno_col].strip()})
    assert len(pheno_ids) > 0, f"no strain IDs read from {args.pheno} column {args.pheno_col}"

    vcf_ids = vcf_sample_ids(args.vcf)
    vcf_norm = {normalize(v): v for v in vcf_ids}

    rows = []
    claims = defaultdict(list)  # vcf_id -> [phenotype_strain_id, ...]

    for pid in pheno_ids:
        vid, tier, score = None, "unmatched", 0.0
        if pid in vcf_ids:
            vid, tier, score = pid, "exact", 1.0
        elif normalize(pid) in vcf_norm:
            vid, tier, score = vcf_norm[normalize(pid)], "normalized", 1.0
        else:
            sfx = suffix_match(pid, vcf_ids)
            if sfx is not None:
                vid, tier, score = sfx, "fuzzy", 0.9
            else:
                fid, fscore = fuzzy_match(pid, vcf_ids)
                if fid is not None:
                    vid, tier, score = fid, "fuzzy", round(fscore, 3)
        rows.append({"phenotype_strain_id": pid, "vcf_sample_id": vid or "", "tier": tier, "match_score": score})
        if vid is not None:
            claims[vid].append(pid)

    collisions = {vid: pids for vid, pids in claims.items() if len(pids) > 1}

    out_dir = args.out_dir
    import os
    os.makedirs(out_dir, exist_ok=True)

    table_path = os.path.join(out_dir, "strain_match_table.csv")
    with open(table_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["phenotype_strain_id", "vcf_sample_id", "tier", "match_score"])
        writer.writeheader()
        writer.writerows(rows)

    review_path = os.path.join(out_dir, "NEEDS_REVIEW.md")
    with open(review_path, "w") as f:
        f.write("# Strain reconciliation — needs human review\n\n")
        if collisions:
            f.write("## COLLISIONS (hard error — must resolve before proceeding)\n\n")
            for vid, pids in collisions.items():
                f.write(f"- VCF sample `{vid}` claimed by {len(pids)} phenotype strains: {pids}\n")
            f.write("\n")
        fuzzy_rows = [r for r in rows if r["tier"] == "fuzzy"]
        unmatched_rows = [r for r in rows if r["tier"] == "unmatched"]
        f.write(f"## Fuzzy matches ({len(fuzzy_rows)}) — verify each by eye\n\n")
        for r in fuzzy_rows:
            f.write(f"- `{r['phenotype_strain_id']}` -> `{r['vcf_sample_id']}` (score={r['match_score']})\n")
        f.write(f"\n## Unmatched ({len(unmatched_rows)}) — no genotype available\n\n")
        for r in unmatched_rows:
            f.write(f"- `{r['phenotype_strain_id']}`\n")

    n_exact = sum(1 for r in rows if r["tier"] == "exact")
    n_norm = sum(1 for r in rows if r["tier"] == "normalized")
    n_fuzzy = sum(1 for r in rows if r["tier"] == "fuzzy")
    n_unmatched = sum(1 for r in rows if r["tier"] == "unmatched")
    print(f"Phenotype strains: {len(pheno_ids)}")
    print(f"  exact: {n_exact}  normalized: {n_norm}  fuzzy: {n_fuzzy}  unmatched: {n_unmatched}")
    print(f"  collisions: {len(collisions)}")
    print(f"Wrote {table_path}")
    print(f"Wrote {review_path}")

    if collisions:
        print("ERROR: collisions found — resolve in NEEDS_REVIEW.md before creating "
              "strain_match_table.reviewed.csv", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
