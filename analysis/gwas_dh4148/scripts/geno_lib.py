"""Shared loaders for the panel genotypes (BIMBAM, haploid 0/1) and the DH4148 annotation."""
import gzip, re
from pathlib import Path
import numpy as np, pandas as pd

R = Path("analysis/gwas_dh4148/results"); GENO = R / "geno"
FAI = "/bigdata/stajichlab/shared/projects/Rhodotorula/PopGen/Rhodotorula_mucilaginosa_DH4148_ref/refgenome/GCA_058775505.1_UCR_RmucDH4148_1.0_genomic.fna.fai"
GFF = "/bigdata/stajichlab/shared/projects/Rhodotorula/PopGen/Rhodotorula_mucilaginosa_DH4148_ref/refgenome/GCA_058775505.1_UCR_RmucDH4148_1.0_genomic.gff"

def load_geno():
    """Returns snp ids, genotype matrix (n_snp x 126, float32, NaN = missing), strain_id array in column order."""
    ids, rows = [], []
    with gzip.open(GENO / "panel.bimbam.gz", "rt") as f:
        for line in f:
            p = line.rstrip("\n").split(",")
            ids.append(p[0]); rows.append([np.nan if x == "NA" else float(x) for x in p[3:]])
    G = np.asarray(rows, dtype=np.float32)
    order = [l.strip() for l in open(GENO / "panel_samples_in_vcf_order.txt")]
    ps = pd.read_csv(GENO / "panel_strains.csv", dtype=str).set_index("popgen_strain").loc[order]
    return np.array(ids), G, ps.strain_id.astype(int).values, ps.sample_name.values

def split_id(ids):
    x = pd.Series(ids).str.split(":", expand=True); return x[0].values, x[1].astype(int).values

def contigs():
    f = pd.read_csv(FAI, sep="\t", header=None, usecols=[0, 1], names=["chr", "len"]); return f[f.len >= 100000].reset_index(drop=True)

def load_genes():
    g = {}; prod = {}
    for l in open(GFF):
        if l.startswith("#"): continue
        c = l.rstrip("\n").split("\t")
        if len(c) < 9: continue
        at = dict(kv.split("=", 1) for kv in c[8].split(";") if "=" in kv)
        if c[2] == "gene": g[at["ID"].replace("gene-", "")] = dict(gene_id=at["ID"].replace("gene-", ""), chr=c[0], start=int(c[3]), end=int(c[4]), strand=c[6], symbol=at.get("gene", ""))
        elif c[2] == "mRNA" and "Parent" in at:
            gid = at["Parent"].replace("gene-", ""); prod.setdefault(gid, re.sub(r"%2C", ",", at.get("product", "")))
    df = pd.DataFrame(g.values()); df["product"] = df.gene_id.map(prod).fillna(""); return df

IMPACT = {"HIGH": 3, "MODERATE": 2, "LOW": 1, "MODIFIER": 0}
def best_ann(ann, alt):
    """Most severe snpEff annotation for the ALT allele of one site: (effect, impact, gene_id, hgvs_p)."""
    best = ("", "", "", ""); bs = -1
    for a in ann.split(","):
        f = a.split("|")
        if len(f) < 11 or f[0] != alt: continue
        s = IMPACT.get(f[2], 0)
        if s > bs: bs = s; best = (f[1], f[2], f[4], f[10])
    return best

def ann_for(snp_ids):
    want = set(snp_ids); out = {}
    with gzip.open(GENO / "panel_sites_ann.tsv.gz", "rt") as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            k = f"{c[0]}:{c[1]}:{c[2]}:{c[3]}"
            if k in want: out[k] = best_ann(c[6], c[3]) if len(c) > 6 else ("", "", "", "")
    return out
