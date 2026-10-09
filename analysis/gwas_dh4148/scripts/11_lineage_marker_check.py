#!/usr/bin/env python3
"""Are the GWAS lead SNPs lineage markers? A SNP is a lineage-level marker when it does not vary inside any lineage (clonal group, <= 2,000 SNP differences).
Adds that flag and the number of lineages carrying the alternate allele to each loci table, and summarises per scan."""
import sys
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, "analysis/gwas_dh4148/scripts")
from geno_lib import *
T = Path("analysis/gwas_dh4148/report/tables")
ids, X, sid, names = load_geno(); ix = {k: i for i, k in enumerate(ids)}
lin = pd.read_csv(R / "lineages.csv").set_index("strain_id").loc[sid, "lineage"].values; ul = np.unique(lin)
def lineage_stats(snp):
    g = X[ix[snp]]; carries = 0; mixed = 0
    for l in ul:
        v = g[lin == l]; v = v[~np.isnan(v)]
        if len(v) == 0: continue
        if v.min() != v.max(): mixed += 1
        if v.mean() >= 0.5: carries += 1
    return carries, mixed
rows = []
for tag in ("_unadjusted", "_runadj", "_runadj_lineage"):
    f = T / f"gwas_loci{tag}.csv"
    if not f.exists(): continue
    d = pd.read_csv(f); st = [lineage_stats(s) for s in d.lead_snp]; d["lineages_carrying_alt"] = [a for a, b in st]; d["lineages_polymorphic_inside"] = [b for a, b in st]; d["lineage_level_marker"] = d.lineages_polymorphic_inside == 0
    d.to_csv(f, index=False)
    b = d[d.bonferroni_sig]
    rows.append(dict(scan=tag.strip("_"), loci_suggestive=len(d), loci_bonferroni=len(b), bonferroni_lineage_markers=int(b.lineage_level_marker.sum()), contigs_with_bonferroni_loci=b.chr.nunique(),
                     median_lineages_carrying_alt=float(b.lineages_carrying_alt.median()) if len(b) else np.nan))
pd.DataFrame(rows).to_csv(T / "gwas_lead_snp_lineage_check.csv", index=False); print(pd.DataFrame(rows).to_string(index=False))
