#!/usr/bin/env python3
"""Greedy near-clone culling by pairwise IBS0 rate, reconstructed from
analysis/ideas/2026-08-15-color-phenotype-space/PROGRESS.md section 6 (N2):
"IBS0<0.005 greedy" -> 173 informative strains from 201. This algorithm was run
ad hoc on $SCRATCH in the original analysis and never saved as a script; this is
the first reusable implementation.

Pairwise IBS0 (proportion of considered SNP pairs that are opposite homozygotes)
comes from plink2 `--make-king-table cols=id,nsnp,ibs0`. Since every strain in
this panel is haploid-encoded (0 het genome-wide, confirmed by check_ploidy.py),
IBS0 here reduces to a genotype-mismatch rate between two haploid-as-homozygous
samples -- exactly the near-clone signal PROGRESS.md describes.

Greedy removal: repeatedly take the pair with the SMALLEST IBS0 (most similar);
if it is below --threshold, drop one member (deterministic tie-break: the
lexicographically later strain ID, since no per-strain "better genotyped"
signal was recorded in the original run to break ties by); repeat until no
remaining pair is below threshold.

Usage:
  pixi run python3 analysis/gwas/scripts/cull_near_clones.py \
      --bfile analysis/gwas/results/gwas/grm_conditioning/rebuilt_kinship/gwas.pruned \
      --plink2 <path-to-plink2> \
      --threshold 0.005 \
      --out analysis/gwas/results/gwas/near_clone_culling/culled_keep.txt
"""
import argparse
import os
import subprocess

import numpy as np
import pandas as pd


def run_king_table(plink2_bin: str, bfile: str, out_prefix: str) -> str:
    cmd = [
        plink2_bin, "--bfile", bfile,
        "--make-king-table", "cols=id,nsnp,ibs0",
        "--out", out_prefix, "--allow-extra-chr",
    ]
    r = subprocess.run(cmd, capture_output=True, text=True)
    kin0 = out_prefix + ".kin0"
    if r.returncode != 0 or not os.path.exists(kin0):
        raise RuntimeError(f"plink2 --make-king-table failed:\n{r.stdout[-2000:]}\n{r.stderr[-2000:]}")
    return kin0


def greedy_cull(pairs: pd.DataFrame, all_ids: list[str], threshold: float):
    """pairs: DataFrame with columns IID1, IID2, IBS0, sorted ascending by IBS0.
    Repeatedly drop one member of the closest below-threshold pair, re-checking
    remaining pairs each time (a dropped strain's other pairs no longer count)."""
    alive = set(all_ids)
    removed_log = []
    p = pairs[["IID1", "IID2", "IBS0"]].copy()
    while True:
        cand = p[p["IID1"].isin(alive) & p["IID2"].isin(alive)]
        if cand.empty:
            break
        row = cand.loc[cand["IBS0"].idxmin()]
        if row["IBS0"] >= threshold:
            break
        a, b = row["IID1"], row["IID2"]
        drop = b if b > a else a
        removed_log.append((a, b, float(row["IBS0"]), drop))
        alive.discard(drop)
    return sorted(alive), removed_log


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bfile", required=True, help="plink bfile root (LD-pruned)")
    ap.add_argument("--plink2", default="plink2", help="path to plink2 binary")
    ap.add_argument("--threshold", type=float, default=0.005)
    ap.add_argument("--restrict-fam", default=None,
                     help="optional: restrict to strain IDs in this .fam (validation mode)")
    ap.add_argument("--out", required=True, help="output path for kept strain list (one ID per line)")
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    work_prefix = os.path.join(os.path.dirname(args.out), "king_table")

    if args.restrict_fam:
        restrict_ids = set(l.split()[1] for l in open(args.restrict_fam) if l.split())
        fam_rows = [l.split() for l in open(args.bfile + ".fam") if l.split()]
        keep_rows = [(fid, iid) for fid, iid, *_ in fam_rows if iid in restrict_ids]
        assert len(keep_rows) == len(restrict_ids), (
            f"restrict-fam has {len(restrict_ids)} IDs but only {len(keep_rows)} found in {args.bfile}.fam"
        )
        subset_prefix = work_prefix + "_subset"
        keep_file = subset_prefix + ".keep.txt"
        with open(keep_file, "w") as f:
            for fid, iid in keep_rows:
                f.write(f"{fid}\t{iid}\n")
        cmd = [args.plink2, "--bfile", args.bfile, "--keep", keep_file,
               "--make-bed", "--out", subset_prefix, "--allow-extra-chr"]
        r = subprocess.run(cmd, capture_output=True, text=True)
        assert r.returncode == 0, f"plink2 --keep failed:\n{r.stdout[-2000:]}"
        bfile_used = subset_prefix
        all_ids = [iid for _, iid in keep_rows]
    else:
        bfile_used = args.bfile
        all_ids = [l.split()[1] for l in open(args.bfile + ".fam") if l.split()]

    kin0 = run_king_table(args.plink2, bfile_used, work_prefix)
    pairs = pd.read_csv(kin0, sep="\t")
    pairs.columns = [c.lstrip("#") for c in pairs.columns]
    assert "IBS0" in pairs.columns, f"IBS0 column missing from {kin0}: {pairs.columns.tolist()}"
    pairs = pairs.sort_values("IBS0").reset_index(drop=True)

    kept, removed_log = greedy_cull(pairs, all_ids, args.threshold)
    with open(args.out, "w") as f:
        f.write("\n".join(kept) + "\n")

    log_path = args.out.replace(".txt", "_removed_log.csv")
    with open(log_path, "w") as f:
        f.write("pair_a,pair_b,ibs0,dropped\n")
        for a, b, d, drop in removed_log:
            f.write(f"{a},{b},{d:.6f},{drop}\n")

    print(f"Input strains: {len(all_ids)}")
    print(f"Removed (near-clone pairs, IBS0 < {args.threshold}): {len(removed_log)}")
    print(f"Kept: {len(kept)}")
    print(f"Wrote {args.out}")
    print(f"Wrote removal log: {log_path}")


if __name__ == "__main__":
    main()
