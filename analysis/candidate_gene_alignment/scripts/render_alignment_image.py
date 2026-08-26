#!/usr/bin/env python3
"""Static color-block alignment PNG for one gene's polymorphic CDS positions across
the 213-strain panel. Per user decision (2026-08-26): static images + CSV only this
round, no interactive viewer (flagged as a follow-up in the analysis doc).

Rows (strains) are sorted either by a phenotype value (pathway genes with no locus)
or by lead-SNP genotype then phenotype (GWAS-locus genes) so the color blocks reveal
whether allele patterns track the trait, not just position in an arbitrary strain
order. Columns get a coordinate ruler (CDS position) per the bioinformatics review.
"""
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BASE_COLORS = {"A": "#4daf4a", "C": "#377eb8", "G": "#ffff33", "T": "#e41a1c", "N": "#bdbdbd"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dna-polymorphic-csv", required=True)
    ap.add_argument("--strain-context-csv", required=True, help="from build_variant_table.py, for sort key")
    ap.add_argument("--sort-by", default="phenotype", choices=["phenotype", "lead_snp"])
    ap.add_argument("--phenotype-col", default="chroma")
    ap.add_argument("--gene-id", required=True)
    ap.add_argument("--out-png", required=True)
    args = ap.parse_args()

    dna = pd.read_csv(args.dna_polymorphic_csv).set_index("strain")
    ctx = pd.read_csv(args.strain_context_csv).set_index("strain")
    ctx = ctx.reindex(dna.index)

    if args.sort_by == "lead_snp" and "lead_snp_genotype" in ctx.columns:
        order_key = ctx["lead_snp_genotype"].fillna("missing")
        sort_cols = [order_key, ctx[args.phenotype_col] if args.phenotype_col in ctx.columns else 0]
        order = pd.DataFrame({"a": order_key, "b": sort_cols[1]}).sort_values(["a", "b"]).index
    else:
        col = args.phenotype_col if args.phenotype_col in ctx.columns else ctx.columns[0]
        order = ctx.sort_values(col).index

    dna = dna.loc[order]
    pos_cols = [c for c in dna.columns]
    positions = [int(c.replace("cds_pos_", "")) for c in pos_cols]

    n_strains, n_pos = dna.shape
    fig_w = max(6, min(24, n_pos * 0.18))
    fig_h = max(4, min(40, n_strains * 0.06))
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    mat = dna.to_numpy()
    color_ids = {b: i for i, b in enumerate(["A", "C", "G", "T", "N"])}
    idx_mat = np.vectorize(lambda b: color_ids.get(b, 4))(mat)
    cmap = matplotlib.colors.ListedColormap([BASE_COLORS[b] for b in ["A", "C", "G", "T", "N"]])
    ax.imshow(idx_mat, aspect="auto", cmap=cmap, vmin=-0.5, vmax=4.5, interpolation="none")

    tick_step = max(1, n_pos // 40)
    ax.set_xticks(range(0, n_pos, tick_step))
    ax.set_xticklabels([positions[i] for i in range(0, n_pos, tick_step)], rotation=90, fontsize=5)
    ax.set_xlabel("CDS position (1-based)")

    ytick_step = max(1, n_strains // 40)
    ax.set_yticks(range(0, n_strains, ytick_step))
    ax.set_yticklabels([order[i] for i in range(0, n_strains, ytick_step)], fontsize=4)
    ax.set_ylabel(f"strain (sorted by {args.sort_by})")
    ax.set_title(f"{args.gene_id}: {n_pos} polymorphic CDS positions, {n_strains} strains")

    legend_handles = [plt.Rectangle((0, 0), 1, 1, color=BASE_COLORS[b]) for b in ["A", "C", "G", "T", "N"]]
    ax.legend(legend_handles, ["A", "C", "G", "T", "N (missing)"], loc="upper right",
              bbox_to_anchor=(1.15, 1.0), fontsize=6, title="allele")

    fig.tight_layout()
    fig.savefig(args.out_png, dpi=200)
    print(f"{args.gene_id}: alignment image -> {args.out_png} ({n_strains} strains x {n_pos} positions)")


if __name__ == "__main__":
    main()
