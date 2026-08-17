#!/usr/bin/env python3
"""Tier D gene map + Tier E credible-set resolution figure.

Two-panel figure:
  [A] Gene-map around the scaffold_10 anchor region (spans the chroma and
      AUC_10 leads ~11 kb apart): -log10 p of FDR-significant SNPs coloured by
      trait across the window, with the lead SNP marked and a schematic gene
      track below spanning the genes hit by overlapping SNPs (pulled from the
      Tier-D annotated CSV, so gene products are real from the annotation).
  [B] Tier E 95% credible-set resolution: for each of the 14 anchors, the
      size of the 95% credible set (x, log) vs lead posterior probability
      (y), point coloured by the rare-driven flag and annotated with the
      nearest/overlapping gene product.
Saves PNG + PDF.  Uses only Tier-D/Tier-E outputs (no per-SNP assoc needed).
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

LOCUS = dict(traits=("chroma", "AUC_10"), chrom="10",
             x0=364_905, x1=416_172)   # ±20 kb around both scaffold_10 leads

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sig", required=True, help="tierD_fdr_snps_annotated.csv.gz")
    ap.add_argument("--credset", required=True, help="tierE_credible_sets.csv")
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()
    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)

    sig = pd.read_csv(args.sig)
    cs = pd.read_csv(args.credset)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.2))

    # ------------------------- Panel A: gene map -------------------------
    ax = axes[0]
    colours = {"chroma": "#1b9e77", "AUC_10": "#d95f02"}
    win = sig[(sig["chr"].astype(str) == LOCUS["chrom"])
              & (sig["ps"].between(LOCUS["x0"], LOCUS["x1"]))]
    for t in LOCUS["traits"]:
        g = win[win["trait"] == t]
        if len(g) == 0:
            continue
        ax.scatter(g["ps"], -np.log10(g["p_wald"]), s=14, alpha=.75,
                   color=colours[t], edgecolors="none",
                   label=t if t in g["trait"].values else None)
        lead = g.loc[g["p_wald"].idxmin()]
        ax.annotate(f"{lead['rs']}\n−log10p={-np.log10(lead['p_wald']):.1f}",
                    (lead["ps"], -np.log10(lead["p_wald"])),
                    xytext=(6, 6), textcoords="offset points", fontsize=7,
                    arrowprops=dict(arrowstyle="-", lw=.6))
    ax.set_xlim(LOCUS["x0"], LOCUS["x1"])
    ax.set_ylim(0, max(-np.log10(win["p_wald"])) * 1.2 + .5)
    ax.set_xlabel("position on scaffold_10 (bp)")
    ax.set_ylabel("-log10 p (Tier-D FDR-significant SNP)")
    ax.set_title(f"(A) Tier D gene map — chromosome {LOCUS['chrom']} anchor region")
    ax.legend(fontsize=8, loc="upper right")

    # schematic gene track: spans from SNP positions overlapping each gene
    genes = (win[win["overlap_gene"].notna()]
             .groupby("overlap_gene")
             .agg(gstart=("ps", "min"), gend=("ps", "max"),
                  product=("product", "first"), trait=("trait", "first"))
             .reset_index())
    track_ymin, track_ymax = -0.09, 0.04
    ylim_lo = min(0, track_ymin - .06)
    ax.set_ylim(ylim_lo, ax.get_ylim()[1])
    for _, r in genes.iterrows():
        c = colours[r["trait"]] if r["trait"] in colours else "0.5"
        ax.plot([r["gstart"], r["gend"]], [track_ymin, track_ymin],
                color=c, lw=4, solid_capstyle="butt")
        ax.plot(r["gstart"], track_ymin, "|", ms=12, color=c)
        ax.plot(r["gend"], track_ymin, "|", ms=12, color=c)
        ax.text((r["gstart"] + r["gend"]) / 2, track_ymin - .03,
                r["overlap_gene"], ha="center", fontsize=6, rotation=45)
        prod = "hypothetical protein" if r["product"] == "hypothetical protein" \
            else f"{r['product']} ({r['overlap_gene']})"
        ax.text((r["gstart"] + r["gend"]) / 2, track_ymin - .06,
                r["product"], ha="center", fontsize=6, rotation=45,
                color=c)
    ax.axhline(0, color="0.6", lw=.8, ls=":")

    # --------------------- Panel B: credible sets ------------------------
    ax = axes[1]
    p95 = cs[cs["credset"] == "95%"].copy()
    rare_c = {True: "#e31a1c", False: "#3182bd"}
    for _, r in p95.iterrows():
        ax.scatter(r["n_snps_credset"], r["lead_pp"], s=75,
                   color=rare_c[r["rare_driven"]], edgecolors="k", lw=.4,
                   alpha=.9)
        prod = ("DBP3 / RNA-dependent ATPase" if r["nearest_gene"] == "OM429_004640"
                else r["product"])
        if r["gene_rel"] == "inside":
            tag = f"{r['label']} [{r['nearest_gene']}]"
        else:
            tag = f"{r['label']}"
        ax.annotate(tag, (r["n_snps_credset"], r["lead_pp"]),
                    xytext=(4, 5), textcoords="offset points", fontsize=6.5)
    ax.set_xscale("log")
    ax.set_xlabel("95% credible-set size (n SNPs)")
    ax.set_ylabel("lead posterior probability (ABF, z-space)")
    ax.set_title("(B) Tier E credible-set resolution (14 anchors, 95% CS)")
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([], [], marker="o", ls="", color=rare_c[True],
                              label="rare-driven (AF<0.05)"),
                       Line2D([], [], marker="o", ls="", color=rare_c[False],
                              label="common variant")], fontsize=8, loc="lower left")

    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(out / f"tierde_gene_finemap.{ext}", dpi=150)
    print(f"wrote tierde_gene_finemap.{{png,pdf}}  "
          f"(A: {len(win)} sig SNPs, {len(genes)} genes in window; "
          f"B: {len(p95)} anchors)")

if __name__ == "__main__":
    main()
