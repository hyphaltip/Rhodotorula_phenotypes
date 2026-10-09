#!/usr/bin/env python3
"""Cross-species homolog comparison table for the candidate genes' significant
coding variants: for each FDR-significant/likely_real missense variant, report
the top diamond blastp hit's identity/product in 10 other Rhodotorula species +
2 Cystobasidium outgroups + S. cerevisiae, and whether the mutated residue
position is conserved (same AA), similar, or diverged in that homolog.

Inputs:
  analysis/candidate_gene_alignment/results/candidate_gene_phenotype_assoc_all10.csv
  analysis/gwas/results/candidate_genes/homology/blast/*.tsv  (diamond outfmt 6 with qseq/sseq)
Output:
  analysis/gwas/results/candidate_genes/homology/homolog_impact_table.csv
  analysis/gwas/results/candidate_genes/homology/homolog_impact_table.md
"""
import re
import glob
import os
import pandas as pd

CGA = "analysis/candidate_gene_alignment/results/candidate_gene_phenotype_assoc_all10.csv"
BLASTDIR = "analysis/gwas/results/candidate_genes/homology/blast"
OUTDIR = "analysis/gwas/results/candidate_genes/homology"

GENE_PRODUCT = {
    "OM429_003333": "phytoene synthase/lycopene cyclase (psy/lcy)",
    "OM429_003336": "phytoene desaturase (pds)",
    "OM429_000065": "sat-locus gene (hypothetical)",
    "OM429_001415": "lab_L GWAS locus gene",
    "OM429_001430": "lab_L GWAS locus gene",
    "OM429_001521": "8-oxoguanine glycosylase (ogg1)",
    "OM429_001533": "cu_dose GWAS locus gene",
    "OM429_002663": "saccharopine dehydrogenase (LYS1)",
    "OM429_003729": "histone chaperone (NAP1)",
    "OM429_005034": "GTPase-activating protein (GYP1)",
}

AA3to1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q",
    "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
    "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
    "Tyr": "Y", "Val": "V",
}

HGVS_RE = re.compile(r"p\.([A-Za-z]{3})(\d+)([A-Za-z]{3})")

COLS = ["qseqid", "sseqid", "pident", "length", "evalue", "bitscore",
        "qstart", "qend", "sstart", "send", "qlen", "slen", "qseq", "sseq", "stitle"]


def load_blast():
    hits = {}
    for f in glob.glob(f"{BLASTDIR}/*.tsv"):
        species = os.path.basename(f).replace(".tsv", "")
        df = pd.read_csv(f, sep="\t", header=None, names=COLS)
        for _, row in df.iterrows():
            gene = row["qseqid"].split("|")[0]
            hits.setdefault(gene, {})[species] = row
    return hits


def aligned_residue(row, query_pos_1based):
    """Map a 1-based query protein position to the aligned subject residue."""
    qstart, qend = int(row["qstart"]), int(row["qend"])
    if query_pos_1based < qstart or query_pos_1based > qend:
        return None, None  # outside the aligned region
    qseq, sseq = row["qseq"], row["sseq"]
    q_ungapped = qstart - 1
    for qc, sc in zip(qseq, sseq):
        if qc != "-":
            q_ungapped += 1
        if q_ungapped == query_pos_1based:
            return qc, sc
    return None, None


def main():
    assoc = pd.read_csv(CGA)
    sig = assoc[(assoc["fdr_sig"]) & (assoc["verdict"] == "likely_real") &
                (assoc["consequence"] == "missense_variant")].copy()
    variants = sig.drop_duplicates(subset=["gene", "pos", "hgvs_p"])[["gene", "pos", "hgvs_p", "hgvs_c"]]

    blast = load_blast()
    species_order = sorted(blast[next(iter(blast))].keys()) if blast else []

    rows = []
    for _, v in variants.iterrows():
        gene, hgvs_p = v["gene"], v["hgvs_p"]
        m = HGVS_RE.match(hgvs_p)
        if not m:
            continue
        ref3, pos, alt3 = m.groups()
        pos = int(pos)
        ref1, alt1 = AA3to1.get(ref3, "?"), AA3to1.get(alt3, "?")

        for species, row in blast.get(gene, {}).items():
            q_res, s_res = aligned_residue(row, pos)
            if q_res is None:
                status = "outside_aligned_region"
            elif q_res != ref1:
                status = f"coordinate_mismatch(query={q_res})"
            elif s_res == "-":
                status = "gap_in_homolog"
            elif s_res == ref1:
                status = "conserved(matches R. mucilaginosa ref)"
            else:
                status = f"diverged(homolog={s_res})"

            rows.append({
                "gene": gene, "product": GENE_PRODUCT.get(gene, ""),
                "hgvs_c": v["hgvs_c"], "hgvs_p": hgvs_p,
                "species": species.replace("_", " "),
                "homolog_id": row["sseqid"], "pident_full_protein": row["pident"],
                "evalue": row["evalue"], "bitscore": row["bitscore"],
                "homolog_product": row["stitle"],
                "residue_status_at_variant": status,
            })

    out = pd.DataFrame(rows)
    out.to_csv(f"{OUTDIR}/homolog_impact_table.csv", index=False)

    # Markdown summary: one row per gene x variant, columns = conservation call per species (compact)
    if not out.empty:
        pivot = out.pivot_table(index=["gene", "hgvs_p"], columns="species",
                                 values="residue_status_at_variant", aggfunc="first")
        pivot = pivot.map(lambda s: {
            "conserved(matches R. mucilaginosa ref)": "conserved",
        }.get(s, s) if isinstance(s, str) else s)
        pivot = pivot.map(lambda s: "diverged" if isinstance(s, str) and s.startswith("diverged") else s)
        pivot = pivot.fillna("no hit")

        def write_markdown_table(df, fh):
            cols = list(df.columns)
            fh.write("| gene | variant | " + " | ".join(cols) + " |\n")
            fh.write("|---|---|" + "---|" * len(cols) + "\n")
            for (gene, hgvs_p), r in df.iterrows():
                fh.write(f"| {gene} | {hgvs_p} | " + " | ".join(str(r[c]) for c in cols) + " |\n")

        with open(f"{OUTDIR}/homolog_impact_table.md", "w") as fh:
            write_markdown_table(pivot, fh)

        # Compact summary: conservation fraction across aligned homologs + S. cerevisiae call + best Rhodotorula homolog product
        summary_rows = []
        for (gene, hgvs_p), sub in out.groupby(["gene", "hgvs_p"]):
            aligned = sub[~sub["residue_status_at_variant"].isin(
                ["outside_aligned_region"]) & ~sub["residue_status_at_variant"].str.startswith("coordinate_mismatch")]
            n_conserved = (aligned["residue_status_at_variant"] == "conserved(matches R. mucilaginosa ref)").sum()
            n_aligned = len(aligned)
            scer = sub[sub["species"] == "Saccharomyces cerevisiae"]
            scer_call = scer["residue_status_at_variant"].iloc[0] if not scer.empty else "no hit"
            paralog = sub[sub["species"] == "Rhodotorula mucilaginosa NRRL Y-2510"]
            cross = sub[sub["species"] != "Rhodotorula mucilaginosa NRRL Y-2510"]
            best = cross.loc[cross["bitscore"].idxmax()] if not cross.empty else None
            summary_rows.append({
                "gene": gene, "product": GENE_PRODUCT.get(gene, ""), "hgvs_p": hgvs_p,
                "conserved_in_N_of_M_aligned_homologs": f"{n_conserved}/{n_aligned}",
                "S_cerevisiae_residue_call": scer_call,
                "same_genome_paralog": paralog["homolog_id"].iloc[0] if not paralog.empty else "none",
                "top_cross_species_homolog": best["species"] if best is not None else "n/a",
                "top_cross_species_pident_%": round(best["pident_full_protein"], 1) if best is not None else None,
                "top_cross_species_product": best["homolog_product"] if best is not None else "n/a",
            })
        pd.DataFrame(summary_rows).to_csv(f"{OUTDIR}/homolog_impact_summary.csv", index=False)
    print(f"wrote {len(out)} rows across {variants.shape[0]} variants x up to {len(species_order)} species")


if __name__ == "__main__":
    main()
