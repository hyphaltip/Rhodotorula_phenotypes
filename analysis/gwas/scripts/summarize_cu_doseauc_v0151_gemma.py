#!/usr/bin/env python3
"""Summarize the cu_doseauc_v0151 GEMMA scan: BH-FDR, top hits, and nearest/
overlapping gene annotation (funannotate GFF3), plus a lookup of each hit's
p-value in every other Tier-A trait's assoc file (replication/novelty check).

Usage: pixi run python3 analysis/gwas/scripts/summarize_cu_doseauc_v0151_gemma.py
"""
from __future__ import annotations

import gzip
import pathlib

import pandas as pd
from statsmodels.stats.multitest import multipletests

REPO = pathlib.Path(__file__).resolve().parents[3]
GWAS = REPO / "analysis/gwas"
ASSOC = GWAS / "results/gwas/tierA_summary/gemma_output/gwas_cu_doseauc_v0151.assoc.txt"
OUTDIR = GWAS / "results/gwas/tierA_summary"
GFF3 = pathlib.Path(
    "/bigdata/stajichlab/shared/projects/Population_Genomics/"
    "Rhodotorula_mucilaginosa_NRRLY2510/genome/"
    "Rhodotorula_mucilaginosa_NRRL_Y-2510.gff3.gz"
)
OTHER_TRAITS = [
    "chroma", "sat", "bright", "clone_mean_area",
    "AUC_0", "AUC_10", "AUC_20", "AUC_30", "AUC_ratio_10",
    "resilience_30", "cu_dose_slope", "IC50_est",
]
FDR_ALPHA = 0.05
FLANK = 2000  # bp; matches this repo's candidate-gene CDS+/-2kb convention


def load_genes(gff3_path: pathlib.Path) -> pd.DataFrame:
    rows = []
    with gzip.open(gff3_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            scaffold, start, end, attrs = parts[0], int(parts[3]), int(parts[4]), parts[8]
            gene_id = next((kv.split("=", 1)[1] for kv in attrs.split(";") if kv.startswith("ID=")), None)
            name = next((kv.split("=", 1)[1] for kv in attrs.split(";") if kv.startswith("Name=")), "")
            rows.append({"scaffold": scaffold, "start": start, "end": end,
                         "gene_id": gene_id, "name": name})
    genes = pd.DataFrame(rows)
    # ANALYSIS_OK[runtime-assert]: developer tripwire on the shared GFF3's
    # expected content, checked once at load time.
    assert len(genes) > 0, f"no gene features parsed from {gff3_path}"
    return genes


def nearest_gene(scaffold: str, pos: int, genes: pd.DataFrame) -> str:
    g = genes[genes.scaffold == scaffold]
    inside = g[(g.start - FLANK <= pos) & (g.end + FLANK >= pos)]
    if len(inside):
        inside = inside.assign(
            dist=lambda d: (d.start - pos).abs().where(pos < d.start, 0)
            .where(pos <= d.end, (pos - d.end).abs())
        )
        best = inside.sort_values("dist").iloc[0]
        # NB: use best["name"], not best.name -- the latter is pandas' reserved
        # Series.name (the row index), not this "name" column.
        label = best["gene_id"] + (f" ({best['name']})" if best["name"] else "")
        inside_flag = "inside" if best.start <= pos <= best.end else f"{int(best.dist)}bp"
        return f"{label} [{inside_flag}]"
    return "no gene within 2kb"


def lookup_other_trait_p(scaffold: str, pos: int, trait: str) -> float | None:
    f = OUTDIR / "gemma_output" / f"gwas_{trait}.assoc.txt"
    if not f.exists():
        # ANALYSIS_OK[optional-input]: not every Tier-A trait necessarily has an
        # assoc file in every run of this scan; None here means "not compared",
        # rendered as an explicit blank cell, distinct from "no hit found".
        return None
    df = pd.read_csv(f, sep="\t", usecols=["chr", "ps", "p_wald"])
    chr_num = int(scaffold.replace("scaffold_", ""))
    # ANALYSIS_OK[sample-filter]: intentional locus-window lookup (this hit's
    # position +/- FLANK), not an incidental drop -- the whole point of the
    # function is to test one specific window in the other trait's scan.
    hit = df[(df.chr == chr_num) & (df.ps.between(pos - FLANK, pos + FLANK))]
    if hit.empty:
        return None
    return float(hit.p_wald.min())


def main() -> None:
    # ANALYSIS_OK[sample-filter]: GEMMA writes NA p_wald for SNPs it could not
    # fit (e.g. singular design at that site); excluded from FDR/top-hit
    # ranking since they carry no test statistic, consistent with the -lmm 4
    # run's own "analyzed SNPs" count logged in the .log file.
    df = pd.read_csv(ASSOC, sep="\t").dropna(subset=["p_wald"])
    df["fdr_q"] = multipletests(df["p_wald"], method="fdr_bh")[1]
    n_fdr = int((df.fdr_q < FDR_ALPHA).sum())
    print(f"[summarize_v0151] {len(df)} SNPs tested, {n_fdr} FDR<{FDR_ALPHA}")

    genes = load_genes(GFF3)
    top = df.sort_values("p_wald").head(10).copy()
    top["scaffold"] = "scaffold_" + top["chr"].astype(str)
    top["nearest_gene"] = [nearest_gene(s, int(p), genes) for s, p in zip(top.scaffold, top.ps)]
    for trait in OTHER_TRAITS:
        top[f"minp_{trait}"] = [
            lookup_other_trait_p(s, int(p), trait) for s, p in zip(top.scaffold, top.ps)
        ]

    out_cols = ["scaffold", "ps", "af", "beta", "p_wald", "fdr_q", "nearest_gene"] + \
               [f"minp_{t}" for t in OTHER_TRAITS]
    out = top[out_cols]
    out_path = OUTDIR / "cu_doseauc_v0151_top_hits_annotated.csv"
    out.to_csv(out_path, index=False)
    print(out.to_string(index=False))
    print(f"\n[summarize_v0151] wrote {out_path}")


if __name__ == "__main__":
    main()
