#!/usr/bin/env python3
"""Compare DuckDB strain_info against the old Copper.Strain_info.csv strain set,
the VCF sample IDs, the reviewed match table and the GWAS .fam. Read-only; writes
only to the --out-dir given. Run inside SLURM (needs the DB snapshot logic)."""
import argparse, pathlib, subprocess, sys
import pandas as pd
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import duckdb_inputs as DI

REPO = DI.REPO
ap = argparse.ArgumentParser(); ap.add_argument("--out-dir", required=True); a = ap.parse_args()
out = pathlib.Path(a.out_dir); out.mkdir(parents=True, exist_ok=True)
con = DI.connect_readonly(); new = DI.strain_table(con)
old = pd.read_csv(subprocess.run(["git","show","data-v1-pre-metal-replace:data/metadata/Copper.Strain_info.csv"],
        cwd=REPO, capture_output=True, text=True, check=True).stdout and
        pathlib.Path("/dev/stdin"), nrows=0) if False else None
import io
old = pd.read_csv(io.StringIO(subprocess.run(["git","show","data-v1-pre-metal-replace:data/metadata/Copper.Strain_info.csv"],
        cwd=REPO, capture_output=True, text=True, check=True).stdout))
vcf = [l.strip() for l in open(REPO/"analysis/gwas/results/duckdb_retool/vcf_sample_ids.txt") if l.strip()]
rev = pd.read_csv(REPO/"analysis/gwas/results/strain_reconciliation/strain_match_table.reviewed.csv")
accepted = set(rev[rev.decision=="accept"].phenotype_strain_id)
fam = [l.split()[1] for l in open(REPO/"analysis/gwas/results/gwas/grm_conditioning/rebuilt_kinship/gwas.pruned.fam") if l.split()]
oldS = set(old["Strain"].dropna().str.strip()); newS = set(new.Strain)
print(f"old Strain_info rows={len(old)} distinct Strain={len(oldS)}; new non-control strains={len(new)} distinct sample_name={len(newS)}")
print(f"old∩new={len(oldS&newS)} old-only={len(oldS-newS)} new-only={len(newS-oldS)}")
print("old-only:", sorted(oldS-newS)[:30]); print("new-only:", sorted(newS-oldS)[:30])
print(f"VCF samples={len(vcf)}; new sample_name exact in VCF={len(newS&set(vcf))}; old Strain exact in VCF={len(oldS&set(vcf))}")
print(f"reviewed accepted={len(accepted)}; accepted ⊂ new sample_names: {len(accepted&newS)}/{len(accepted)}; missing={sorted(accepted-newS)}")
print(f"fam strains={len(fam)}; in new sample_name={len(set(fam)&newS)}/{len(fam)}; missing={sorted(set(fam)-newS)}")
rm = new[new.Strain.isin(fam)]
print(f"fam rows in strain_info (strain_id level)={len(rm)}; duplicated sample names among them={rm.Strain.duplicated(keep=False).sum()}")
print("species of fam strains:\n", rm.Species.value_counts().to_string())
print("fam strains: metals_tested counts:\n", rm.metals_tested.value_counts().to_string())
# species check old vs new
oldsp = old.dropna(subset=["Strain"]).drop_duplicates("Strain").set_index("Strain")["Species"]
j = new.set_index("Strain")[["Species"]].join(oldsp.rename("old_species"), how="inner")
j = j[~j.index.duplicated()]
print(f"species agree old vs new for common strains: {(j.Species==j.old_species).sum()}/{len(j)}")
print(j[j.Species!=j.old_species].to_string())
new.to_csv(out/"strain_info_from_duckdb.csv", index=False)
