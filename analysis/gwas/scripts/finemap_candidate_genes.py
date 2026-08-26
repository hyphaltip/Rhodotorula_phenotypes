#!/usr/bin/env python3
"""Fine-map the 8 validated `likely_real` loci (D-17/D-18, GWAS.md S13/S14) toward
candidate genes.

Ported from analysis/ideas/2026-08-15-color-phenotype-space/scripts/{finemap_credible_sets.py,
annotate_gwas_loci.py} (that folder's own Tier D/E, kept as read-only reference; the
methodology -- Wakefield ABF in z-space, candidate filter p<1e-3, 90/95/99% credible
sets by posterior probability -- is unchanged and cited there, not re-derived) and
adapted for this port's specific, already-validated locus set rather than a genome-wide
FDR-significant SNP list:

  1. Wakefield ABF credible sets per locus, same math as the original (z-space,
     prior SD=0.2 on the causal NCP, candidate filter p<1e-3).
  2. Nearest-gene annotation for the lead SNP (as before), PLUS a fuller gene-window
     listing: every gene whose span overlaps the +/-100kb LD block already
     characterized in check_mas_gates.py's LD-decay check (not just the nearest one) --
     since those blocks contain hundreds of SNPs in LD, a single "nearest gene" call is
     too narrow a candidate list for loci with this much local LD.

Reuses the SAME genome annotation index (gene_index.json, 6,799 genes) the original run
built -- this is fixed genome annotation, not project-specific data, so it is referenced
from the ideas folder rather than duplicated.

Output: results/gwas/tierE/candidate_genes_credible_sets.csv (per credible-set-level row,
matching the original tierE_credible_sets.csv shape) and
results/gwas/tierE/candidate_genes_window.csv (every gene in each locus's +/-100kb window,
one row per gene, for a broader candidate list).
"""
import argparse
import json
import math
import os
import sys

import pandas as pd

WINDOW_KB = 100  # matches check_mas_gates.py's LD-decay window
PRIOR_SD = 0.2
W = PRIOR_SD ** 2


def z_abf(beta, se):
    if se is None or se <= 0 or beta is None or beta != beta:
        return -math.inf
    z = beta / se
    r = W / (1.0 + W)
    return 0.5 * math.log(1.0 - r) + 0.5 * r * z * z


def load_genes(path):
    with open(path) as f:
        d = json.load(f)
    by_scaf = {}
    for v in d.values():
        by_scaf.setdefault(v["scaf"], []).append(v)
    for k in by_scaf:
        by_scaf[k].sort(key=lambda x: x["start"])
    return d, by_scaf


def nearest_gene(by_scaf, scaf, pos):
    for g in by_scaf.get(scaf, []):
        if g["start"] <= pos <= g["end"]:
            return g["gene"], 0, "inside"
    hits = [g for g in by_scaf.get(scaf, []) if g["start"] > pos]
    if hits:
        g = min(hits, key=lambda x: x["start"] - pos)
        return g["gene"], g["start"] - pos, "near"
    if by_scaf.get(scaf):
        g = max(by_scaf[scaf], key=lambda x: x["end"])
        return g["gene"], pos - g["end"], "near"
    return None, None, None


def genes_in_window(by_scaf, scaf, lo, hi):
    return [g for g in by_scaf.get(scaf, []) if g["end"] >= lo and g["start"] <= hi]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--loci-csv", required=True, help="population_vs_locus.csv, filtered to likely_real")
    ap.add_argument("--assoc-dir", required=True, help="dir with {panel}_{trait}_assoc.csv.gz")
    ap.add_argument("--gene-index", required=True)
    ap.add_argument("--window-kb", type=int, default=WINDOW_KB)
    ap.add_argument("--cand-p", type=float, default=1e-3)
    ap.add_argument("--out-credsets", required=True)
    ap.add_argument("--out-window-genes", required=True)
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.out_credsets), exist_ok=True)
    genes, by_scaf = load_genes(args.gene_index)
    print(f"Loaded {len(genes)} genes across {len(by_scaf)} scaffolds", file=sys.stderr)

    loci = pd.read_csv(args.loci_csv)
    loci = loci[loci["verdict"] == "likely_real"] if "verdict" in loci.columns else loci
    window = args.window_kb * 1000

    credset_rows = []
    window_gene_rows = []
    assoc_cache = {}

    for _, locus in loci.iterrows():
        trait, panel, rs = locus["trait"], locus["panel"], locus["top_snp"]
        scaf, pos = rs.split(":")[0], int(rs.split(":")[1])
        label = f"{trait}_{panel}_{scaf}_{pos}"

        key = f"{panel}_{trait}_assoc.csv.gz"
        if key not in assoc_cache:
            p = os.path.join(args.assoc_dir, key)
            df = pd.read_csv(p)
            df["scaf"] = "scaffold_" + df["chr"].astype(str)
            assoc_cache[key] = df
        df = assoc_cache[key]

        win = df[(df["scaf"] == scaf) & (df["ps"] >= pos - window) & (df["ps"] <= pos + window)].copy()
        n_window = len(win)
        if n_window == 0:
            print(f"  SKIP {label}: no SNPs in window", file=sys.stderr)
            continue

        # --- Gene-window listing: every gene overlapping the LD block ---
        win_genes = genes_in_window(by_scaf, scaf, pos - window, pos + window)
        for g in win_genes:
            dist_to_lead = 0 if g["start"] <= pos <= g["end"] else min(abs(pos - g["start"]), abs(pos - g["end"]))
            window_gene_rows.append({
                "label": label, "trait": trait, "panel": panel, "top_snp": rs,
                "gene": g["gene"], "gene_start": g["start"], "gene_end": g["end"],
                "dist_to_lead_snp": dist_to_lead, "overlaps_lead_snp": bool(g["start"] <= pos <= g["end"]),
                "product": g.get("product"), "go": ";".join(g.get("go", [])),
                "interpro": ";".join(g.get("interpro", [])), "pfam": ";".join(g.get("pfam", [])),
            })
        print(f"  {label}: {len(win_genes)} genes within +/-{args.window_kb}kb", file=sys.stderr)

        # --- Wakefield ABF credible sets (same method as the original Tier E) ---
        cand = win[win["p_wald"] < args.cand_p].copy()
        if cand.empty:
            print(f"  SKIP {label} credible set: no candidates with p<{args.cand_p}", file=sys.stderr)
            continue
        cand["logABF"] = [z_abf(b, s) for b, s in zip(cand["beta"], cand["se"])]
        cand = cand.sort_values("logABF", ascending=False).reset_index(drop=True)
        finite = cand["logABF"][cand["logABF"] != -math.inf]
        if finite.empty:
            print(f"  SKIP {label} credible set: no finite ABF", file=sys.stderr)
            continue
        lmax = float(cand["logABF"].max())
        lse = lmax if len(finite) == 1 else lmax + math.log(sum(math.exp(x - lmax) for x in finite))
        cand["pp"] = [math.exp(x - lse) if x != -math.inf else 0.0 for x in cand["logABF"]]
        cand["cump"] = cand["pp"].cumsum().round(8)

        gid, dist, rel = nearest_gene(by_scaf, scaf, pos)
        g = genes.get(gid, {}) if gid else {}
        for level, thresh in (("90%", 0.90), ("95%", 0.95), ("99%", 0.99)):
            n_cs = max(1, int((cand["cump"] <= thresh).sum()) + (0 if cand.iloc[0]["cump"] <= thresh else 1))
            n_cs = min(n_cs, len(cand))
            cs_rows = cand.head(n_cs)
            lead_cs = cs_rows.iloc[0]
            low95, high95 = lead_cs["beta"] - 1.96 * lead_cs["se"], lead_cs["beta"] + 1.96 * lead_cs["se"]
            credset_rows.append({
                "label": label, "trait": trait, "panel": panel, "scaffold": scaf, "anchor_pos": pos,
                "n_snps_window": n_window, "n_cand": len(cand),
                "credset": level, "n_snps_credset": n_cs,
                "lead_rs": cand.iloc[0]["rs"], "lead_pos": int(cand.iloc[0]["ps"]),
                "lead_pp": round(float(cand.iloc[0]["pp"]), 4), "lead_af": float(cand.iloc[0]["af"]),
                "lead_beta": float(cand.iloc[0]["beta"]), "lead_se": float(cand.iloc[0]["se"]),
                "beta_95ci_lo": round(float(low95), 4), "beta_95ci_hi": round(float(high95), 4),
                "lead_p_wald": float(cand.iloc[0]["p_wald"]),
                "cs_max_beta": float(cs_rows["beta"].abs().max()),
                "cs_snp_min_pp": float(lead_cs["pp"]),
                "nearest_gene_to_anchor": gid, "gene_dist": dist, "gene_rel": rel,
                "product": g.get("product"), "go": ";".join(g.get("go", [])),
                "rare_driven": bool(cand.iloc[0]["af"] < 0.02),
            })
        print(f"  {label}: 99% CS n={n_cs}, lead {cand.iloc[0]['rs']} pp={float(cand.iloc[0]['pp']):.3f} "
              f"af={float(cand.iloc[0]['af']):.3f} nearest_gene={gid} ({g.get('product')})", file=sys.stderr)

    pd.DataFrame(credset_rows).to_csv(args.out_credsets, index=False)
    pd.DataFrame(window_gene_rows).to_csv(args.out_window_genes, index=False)
    print(f"\nWrote {args.out_credsets}", file=sys.stderr)
    print(f"Wrote {args.out_window_genes}", file=sys.stderr)


if __name__ == "__main__":
    main()
