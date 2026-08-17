#!/usr/bin/env python3
"""Tier D: annotate validated FDR GWAS loci and BSLMM clusters to genes.

Inputs (per trait, in analysis results dir):
  results/gwas/fdr/<trait>_fdr05_sig.txt   -> FDR-significant SNPs (GEMMA format)
  results/gwas/tierC_summary/tierc_bslmm_summary.csv -> BSLMM top-PIP clusters

Outputs (results/gwas/tierD/):
  tierD_fdr_snps_annotated.csv     every FDR SNP with gene overlap/nearest
  tierD_independent_loci.csv       positional clump (250kb) lead SNP per locus
  tierD_bslmm_loci_annotated.csv   BSLMM top-PIP SNPs with gene mapping
"""
import argparse, csv, gzip, json, os, re, sys
from collections import defaultdict

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, '..', 'results', 'gwas'))
TIERD = os.path.join(RES, 'tierD')
FDR = os.path.join(RES, 'fdr')
TCRC = os.path.join(RES, 'tierC_summary', 'tierc_bslmm_summary.csv')

FIELDS = ['chr', 'rs', 'ps', 'n_miss', 'allele1', 'allele0', 'af', 'beta', 'se',
          'logl_H1', 'l_remle', 'l_mle', 'p_wald', 'p_lrt', 'p_score', 'sig']

CLUMP_KB = 250


def load_gene_index(path):
    with open(path) as f:
        d = json.load(f)
    genes = {}
    for gid, v in d.items():
        genes[gid] = v
    by_scaf = defaultdict(list)
    for v in genes.values():
        by_scaf[v['scaf']].append(v)
    for k in by_scaf:
        by_scaf[k].sort(key=lambda x: (x['start'], x['end']))
    return genes, by_scaf


def nearest_gene(by_scaf, scaf, pos):
    """Return (gene_id, distance, relation) for position pos on scaf."""
    hits = by_scaf.get(scaf, [])
    if not hits:
        return None, None, None
    best = None
    bestd = None
    for g in hits:
        if g['start'] <= pos <= g['end']:
            return g['gene'], 0, 'inside'
        d = min(abs(pos - g['start']), abs(pos - g['end']))
        if bestd is None or d < bestd:
            best, bestd = g['gene'], d
    return best, bestd, 'near'


def annotate_snp(genes, by_scaf, row):
    scaf = row['chr']
    if isinstance(scaf, (int, float)) or str(scaf).isdigit():
        scaf = f'scaffold_{int(scaf)}'
    pos = row['ps']
    gid, dist, rel = nearest_gene(by_scaf, scaf, pos)
    g = genes.get(gid, {}) if gid else {}
    return {
        'overlap_gene': gid if rel == 'inside' else None,
        'nearest_gene': gid,
        'gene_dist': dist,
        'gene_rel': rel,
        'product': g.get('product'),
        'go': ';'.join(g.get('go', [])),
        'interpro': ';'.join(g.get('interpro', [])),
        'pfam': ';'.join(g.get('pfam', [])),
    }


def read_fdr(path):
    if not os.path.exists(path):
        return None
    return pd.read_csv(path, sep='\t')


def clump(df, kb=CLUMP_KB):
    """Positional clumping: group by chr, walk sorted SNPs, lead = most significant
    within 250kb that is not already assigned to an older lead."""
    df = df.sort_values(['chr', 'ps']).reset_index(drop=True)
    best = df.sort_values('p_wald').drop_duplicates(subset='chr')
    best = best.reset_index(drop=True)
    lead_id = [False] * len(df)
    for i, r in best.iterrows():
        mask = (df['chr'] == r['chr']) & (abs(df['ps'] - r['ps']) <= kb * 1000)
        idx = df[mask].index
        for j in idx:
            if not lead_id[j]:
                lead_id[j] = True
    df['is_lead'] = lead_id
    return df[df['is_lead']].reset_index(drop=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--index', default=os.path.join(TIERD, 'scripts', 'gene_index.json'))
    ap.add_argument('--fdr', default=FDR)
    ap.add_argument('--tierc', default=TCRC)
    ap.add_argument('--outdir', default=TIERD)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    genes, by_scaf = load_gene_index(args.index)
    print(f"Loaded {len(genes)} genes across {len(by_scaf)} scaffolds", file=sys.stderr)

    out_rows = []
    clumped = []
    for f in sorted(os.listdir(args.fdr)):
        if not f.endswith('_fdr05_sig.txt'):
            continue
        trait = f.replace('_fdr05_sig.txt', '')
        p = os.path.join(args.fdr, f)
        df = read_fdr(p)
        if df is None or df.empty:
            continue
        if 'FDR05_sig' in df.columns:
            df = df[df['FDR05_sig'].astype(str).str.lower() == 'true']
        n = len(df)
        ann = pd.concat([df.reset_index(drop=True),
                         df.apply(lambda r: pd.Series(annotate_snp(genes, by_scaf, r)), axis=1)],
                        axis=1)
        ann.insert(0, 'trait', trait)
        out_rows.append(ann)
        lead = clump(df)
        lead = pd.concat([lead.reset_index(drop=True),
                          lead.apply(lambda r: pd.Series(annotate_snp(genes, by_scaf, r)), axis=1)],
                         axis=1)
        lead.insert(0, 'trait', trait)
        lead.insert(1, 'n_sig', n)
        clumped.append(lead)
        print(f"  {trait}: {n} FDR sig -> {len(lead)} independent loci", file=sys.stderr)

    if out_rows:
        full = pd.concat(out_rows, ignore_index=True)
        full.to_csv(os.path.join(args.outdir, 'tierD_fdr_snps_annotated.csv'), index=False)
    if clumped:
        indep = pd.concat(clumped, ignore_index=True)
        indep.to_csv(os.path.join(args.outdir, 'tierD_independent_loci.csv'), index=False)

    # BSLMM top-PIP clusters
    if os.path.exists(args.tierc):
        bs_rows = []
        with open(args.tierc) as f:
            header = f.readline().rstrip('\n').split(',')
            for line in f:
                fields = csv.reader([line.rstrip('\n')]).__next__()
                d = dict(zip(header, fields))
                trait = d.get('trait', '')
                pipes = re.findall(r'([\w\d_]+):([\d/]+)\s*PIP\s*([\d.]+)', str(d.get('top_PIP_loci_pip1', '')))
                for scaf, poss, pip in pipes:
                    for pstr in poss.split('/'):
                        if pstr:
                            row = {'trait': trait, 'scaffold': scaf, 'pos': int(pstr), 'pip': float(pip)}
                            ann = annotate_snp(genes, by_scaf, {'chr': scaf, 'ps': int(pstr)})
                            row.update(ann)
                            bs_rows.append(row)
        if bs_rows:
            bsdf = pd.DataFrame(bs_rows)
            bsdf.to_csv(os.path.join(args.outdir, 'tierD_bslmm_loci_annotated.csv'), index=False)
            print(f"BSLMM: annotated {len(bs_rows)} top-PIP loci", file=sys.stderr)

    print(f"Wrote outputs to {args.outdir}", file=sys.stderr)


if __name__ == '__main__':
    main()
