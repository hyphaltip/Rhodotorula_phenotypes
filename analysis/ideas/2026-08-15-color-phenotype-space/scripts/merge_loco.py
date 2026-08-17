#!/usr/bin/env python3
"""Merge LOCO per-scaffold GEMMA outputs and compare against Tier-A (full-kinship)
single-SNP scans for the sensitivity analysis.

Produces:
  loco_merged_<prefix>.csv   -- one row per (trait, chr): #snps, GEMMA LOCO lambda_GC,
                                top hit (rs, p, beta), same-scaffold Tier-A comparison
  loco_summary_<prefix>.txt  -- human-readable table

Usage:
  merge_loco.py --loco-dir DIR --tierA-dir DIR --out DIR [--prefix gwas|gwasc]
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

PREFIX_SCAFF = {  # LOCO files: loco_<B>_<trait>_<chr>.assoc.txt
    "ngwas": "gwas", "ngwasc": "gwasc",
}

def gc_lambda(p):
    p = np.asarray(p, dtype=float)
    p = p[np.isfinite(p) & (p > 0)]
    if len(p) < 3:
        return np.nan
    # genomic inflation = median of chi2_1 / theoretical median (0.4549)
    chisq = stats.chi2.isf(p, 1)
    return float(np.median(chisq) / stats.chi2.isf(0.5, 1))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--loco-dir", required=True)
    ap.add_argument("--tierA-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--prefix", default="gwas")
    args = ap.parse_args()
    loco = Path(args.loco_dir); tA = Path(args.tierA_dir); out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    # Which LOCO bfile/outroot matches this prefix's Tier-A set file
    bfile, outroot = (("ngwas", "ngwas") if args.prefix == "gwas"
                      else ("ngwasc", "ngwasc"))
    import re as _re
    traits = sorted({_re.sub(r"_\d+\.assoc\.txt$", "", p.name[len(f"loco_{outroot}_"):])
                     for p in loco.glob(f"loco_{outroot}_*.assoc.txt")})
    print("traits:", traits)

    rows = []
    for t in traits:
        ta = pd.read_csv(tA / f"{args.prefix}_{t}_assoc.csv.gz")
        # Tier-A per-chromosome top + lambda
        ta_c = {c: g for c, g in ta.groupby("chr")}
        for c in range(1, 21):
            fp = loco / f"loco_{outroot}_{t}_{c}.assoc.txt"
            if not fp.exists():
                continue
            lo = pd.read_csv(fp, sep="\t")
            if lo.empty or len(lo) < 3:
                print(f"  {t} chr{c}: empty/short ({len(lo)} rows); skip")
                continue
            lam_lo = gc_lambda(lo["p_wald"].to_numpy() if "p_wald" in lo.columns
                               else lo["p_score"].to_numpy())
            hi_lo = lo.loc[lo["p_wald"].idxmin()]
            # Tier-A on that chromosome (same rs names => join by rs)
            tc = ta_c.get(c)
            if tc is None or len(tc) == 0:
                lam_ta, hi_ta, top_rs_ta, p_ta = np.nan, None, "NA", np.nan
            else:
                lam_ta = gc_lambda(tc["p_wald"].to_numpy())
                h = tc.loc[tc["p_wald"].idxmin()]
                lam_ta, hi_ta, top_rs_ta, p_ta = lam_ta, h, h["rs"], h["p_wald"]
            rows.append(dict(
                trait=t, chr=c, n_snps_LOCO=len(lo),
                lambda_LOCO=round(lam_lo, 3), lambda_TierA_chr=round(lam_ta, 3) if lam_ta == lam_ta else np.nan,
                top_rs_LOCO=hi_lo["rs"], top_p_LOCO=float(hi_lo["p_wald"]), top_beta_LOCO=float(hi_lo["beta"]),
                top_rs_TierA=top_rs_ta, top_p_TierA=float(p_ta) if p_ta == p_ta else np.nan,
            ))
    res = pd.DataFrame(rows)
    res.to_csv(out / f"loco_merged_{args.prefix}.csv", index=False)
    print(f"wrote {len(res)} rows -> loco_merged_{args.prefix}.csv")
    if len(res):
        g = res.groupby("trait")["lambda_LOCO"].agg(["size", "median", "min", "max"])
        print("\nLOCO lambda by trait (across 20 chr):")
        print(g.to_string())
        med_lo = res["lambda_LOCO"].median(); med_ta = res["lambda_TierA_chr"].median()
        print(f"\nmedian lambda: LOCO={med_lo:.3f}  TierA-per-chr={med_ta:.3f}")

if __name__ == "__main__":
    main()
