#!/usr/bin/env python3
"""Manhattan-style significance plots for the candidate-gene fine-mapping results.

Two figures:
  1. candidate_gene_association_manhattan.png -- per-gene, per-variant -log10(meta_p)
     across all 6 color traits (SNP coding variants + indel genotypes), the primary
     "does this candidate gene's own coding variation associate with color" plot.
  2. genome_wide_manhattan_with_candidates.png -- true genome-wide GEMMA Manhattan
     plots (lab_a, chroma, bright) with the 10 candidate-gene loci highlighted, showing
     these genes sit at (or near) genome-wide peaks rather than being an arbitrary pick.

Inputs (already produced by the candidate-gene-alignment pipeline / GWAS Tier A):
  analysis/candidate_gene_alignment/results/candidate_gene_phenotype_assoc_all10.csv
  analysis/candidate_gene_alignment/results/candidate_indel_phenotype_assoc_all10.csv
  analysis/gwas/results/gwas/tierA_summary/assoc_csv/gwas_{lab_a,chroma,bright}_assoc.csv.gz
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

CGA = "analysis/candidate_gene_alignment/results"
GWAS = "analysis/gwas/results/gwas/tierA_summary/assoc_csv"
OUTDIR = "analysis/gwas/figures"

GENE_COORDS = {
    "OM429_003333": ("scaffold_7", 168525, 170261),
    "OM429_003336": ("scaffold_7", 173916, 175999),
    "OM429_000065": ("scaffold_1", 207031, 212775),
    "OM429_001415": ("scaffold_3", 227541, 229892),
    "OM429_001430": ("scaffold_3", 264409, 266988),
    "OM429_001521": ("scaffold_3", 502667, 504489),
    "OM429_001533": ("scaffold_3", 544610, 546394),
    "OM429_002663": ("scaffold_5", 882726, 884304),
    "OM429_003729": ("scaffold_8", 38699, 41014),
    "OM429_005034": ("scaffold_11", 608768, 612742),
}
GENE_LABEL = {
    "OM429_000065": "OM429_000065\n(sat-locus)",
    "OM429_005034": "OM429_005034\n(GYP1)",
    "OM429_001521": "OM429_001521\n(ogg1)",
    "OM429_002663": "OM429_002663\n(LYS1)",
    "OM429_003729": "OM429_003729\n(NAP1)",
    "OM429_001533": "OM429_001533\n(cu_dose)",
    "OM429_001415": "OM429_001415",
    "OM429_001430": "OM429_001430",
    "OM429_003333": "OM429_003333\n(psy/lcy)",
    "OM429_003336": "OM429_003336\n(pds)",
}
TRAIT_COLORS = {
    "lab_L": "#4C72B0", "lab_a": "#DD8452", "lab_b": "#55A868",
    "chroma": "#C44E52", "sat": "#8172B2", "bright": "#937860",
}
GENE_ORDER = list(GENE_COORDS.keys())


def load_candidate_tests():
    snp = pd.read_csv(f"{CGA}/candidate_gene_phenotype_assoc_all10.csv")
    snp["variant_class"] = "SNP"
    indel = pd.read_csv(f"{CGA}/candidate_indel_phenotype_assoc_all10.csv")
    indel["variant_class"] = "indel"
    cols = ["gene", "pos", "trait", "meta_p", "meta_p_fdr", "fdr_sig", "verdict", "variant_class"]
    df = pd.concat([snp[cols], indel[cols]], ignore_index=True)
    df["neg_log10_p"] = -np.log10(df["meta_p"].clip(lower=1e-300))
    return df


def plot_candidate_gene_manhattan(df):
    fig, axes = plt.subplots(2, 5, figsize=(22, 8), sharey=True)
    axes = axes.flatten()
    fdr_thresh = df.loc[df["fdr_sig"], "neg_log10_p"].min() if df["fdr_sig"].any() else None

    for ax, gene in zip(axes, GENE_ORDER):
        sub = df[df["gene"] == gene]
        if sub.empty:
            ax.set_title(GENE_LABEL.get(gene, gene), fontsize=9)
            ax.axis("off")
            continue
        for trait, color in TRAIT_COLORS.items():
            tsub = sub[sub["trait"] == trait]
            if tsub.empty:
                continue
            sig = tsub[tsub["fdr_sig"] & (tsub["verdict"] == "likely_real")]
            nonsig = tsub[~(tsub["fdr_sig"] & (tsub["verdict"] == "likely_real"))]
            ax.scatter(nonsig["pos"], nonsig["neg_log10_p"], s=14, color=color, alpha=0.25,
                       marker="o", linewidths=0)
            ax.scatter(sig["pos"], sig["neg_log10_p"], s=45, color=color, alpha=0.95,
                       marker="D", edgecolors="black", linewidths=0.4, label=trait)
        if fdr_thresh is not None:
            ax.axhline(fdr_thresh, color="grey", linestyle="--", linewidth=0.8)
        ax.set_title(GENE_LABEL.get(gene, gene), fontsize=9)
        ax.tick_params(axis="x", labelrotation=45, labelsize=7)
        ax.set_xlabel("")

    for ax in axes[:5]:
        pass
    fig.text(0.5, 0.02, "position (scaffold coordinate)", ha="center", fontsize=11)
    fig.text(0.08, 0.5, r"$-\log_{10}(\mathrm{meta\ }p)$", va="center", rotation="vertical", fontsize=11)

    handles = [Line2D([0], [0], marker="o", color=c, linestyle="", markersize=6, label=t)
               for t, c in TRAIT_COLORS.items()]
    handles.append(Line2D([0], [0], marker="D", color="black", linestyle="", markersize=7,
                           markerfacecolor="white", label="FDR-sig & likely_real"))
    fig.legend(handles=handles, loc="upper center", ncol=8, bbox_to_anchor=(0.5, 1.02), fontsize=8, frameon=False)
    fig.suptitle(
        "Candidate-gene phenotype association: SNP coding variants (circles) + indel genotypes,\n"
        "faint = not FDR-significant, diamond = FDR-sig & likely_real (population-replicated)",
        fontsize=11, y=1.09)
    fig.tight_layout(rect=[0.08, 0.03, 1, 0.92])
    fig.savefig(f"{OUTDIR}/candidate_gene_association_manhattan.png", dpi=180, bbox_inches="tight")
    fig.savefig(f"{OUTDIR}/candidate_gene_association_manhattan.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_genome_wide(traits=("lab_a", "chroma", "bright")):
    scaffold_order = [f"scaffold_{i}" for i in range(1, 30)]
    fig, axes = plt.subplots(len(traits), 1, figsize=(16, 3.2 * len(traits)), sharex=True)
    if len(traits) == 1:
        axes = [axes]

    offsets = {}
    cum = 0
    max_pos = {}
    for trait in traits:
        d = pd.read_csv(f"{GWAS}/gwas_{trait}_assoc.csv.gz")
        for chrom, sub in d.groupby("chr"):
            key = f"scaffold_{chrom}"
            max_pos[key] = max(max_pos.get(key, 0), sub["ps"].max())
    cum = 0
    for s in scaffold_order:
        if s in max_pos:
            offsets[s] = cum
            cum += max_pos[s] + 5000

    for ax, trait in zip(axes, traits):
        d = pd.read_csv(f"{GWAS}/gwas_{trait}_assoc.csv.gz")
        d["scaffold"] = "scaffold_" + d["chr"].astype(str)
        d["gpos"] = d["scaffold"].map(offsets).astype(float) + d["ps"]
        d["neg_log10_p"] = -np.log10(d["p_wald"].clip(lower=1e-300))
        colors = ["#4C72B0" if i % 2 == 0 else "#8C8C8C" for i in range(len(scaffold_order))]
        cmap = dict(zip(scaffold_order, colors))
        d["color"] = d["scaffold"].map(cmap)
        ax.scatter(d["gpos"], d["neg_log10_p"], s=2, c=d["color"].to_numpy(), alpha=0.5, linewidths=0)

        for gene, (scaf, start, end) in GENE_COORDS.items():
            if scaf not in offsets:
                continue
            gp = offsets[scaf] + (start + end) / 2
            ax.axvline(gp, color="red", linestyle=":", linewidth=0.9, alpha=0.8)
            ax.annotate(gene.replace("OM429_", ""), xy=(gp, ax.get_ylim()[1]), xytext=(0, 2),
                        textcoords="offset points", rotation=90, fontsize=6, color="red", ha="center", va="bottom")

        ax.set_ylabel(f"{trait}\n" + r"$-\log_{10}(p)$", fontsize=9)
        ax.set_xlim(0, cum)

    tick_pos = [offsets[s] + max_pos[s] / 2 for s in scaffold_order if s in offsets]
    tick_lab = [s.replace("scaffold_", "") for s in scaffold_order if s in offsets]
    axes[-1].set_xticks(tick_pos)
    axes[-1].set_xticklabels(tick_lab, fontsize=7)
    axes[-1].set_xlabel("scaffold")
    fig.suptitle("Genome-wide GEMMA association (Tier A) with the 10 candidate-gene loci highlighted (red)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(f"{OUTDIR}/genome_wide_manhattan_with_candidates.png", dpi=170, bbox_inches="tight")
    fig.savefig(f"{OUTDIR}/genome_wide_manhattan_with_candidates.pdf", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    df = load_candidate_tests()
    plot_candidate_gene_manhattan(df)
    plot_genome_wide()
    print("wrote figures to", OUTDIR)
