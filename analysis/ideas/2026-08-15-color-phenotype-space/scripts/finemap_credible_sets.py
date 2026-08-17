#!/usr/bin/env python3
"""Tier E: fine-map validated GWAS loci with Wakefield ABF credible sets + effect sizes.

For each validated anchor locus (from Tier A lead + Tier C BSLMM top-PIP):
  1. Pull all SNPs in the flanking window from the Tier A GEMMA assoc CSV.
  2. Compute approximate Bayes factor (Wakefield 2009) from beta, se, af.
  3. Build 90/95/99% credible sets by posterior probability.
  4. Report per-locus effect size (beta, SE, 95% CI), allele freq, annotation.
  5. Flag rare-variant-driven loci (lead AF < 0.02) as fragile.

Output: results/gwas/tierE/tierE_credible_sets.csv
"""
import argparse, json, math, os, re, sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, '..', 'results', 'gwas'))
TIERE = os.path.join(RES, 'tierE')
ASSOC = os.path.join(RES, 'tierA_summary')
INDEX = os.path.join(RES, 'tierD', 'scripts', 'gene_index.json')

# (trait, scaf, pos, label)
ANCHORS = [
    ('chroma',         'scaffold_10', 384905, 'chroma_lead_10'),
    ('chroma',         'scaffold_3',  26959,  'chroma_bslmm_3a'),
    ('chroma',         'scaffold_3',  137161, 'chroma_bslmm_3b'),
    ('chroma',         'scaffold_3',  132903, 'chroma_bslmm_3c'),
    ('AUC_10',         'scaffold_10', 396172, 'auc10_lead_10'),
    ('AUC_10',         'scaffold_5',  533603, 'auc10_bslmm_5'),
    ('AUC_10',         'scaffold_2',  939277, 'auc10_bslmm_2'),
    ('resilience_30',  'scaffold_13', 810026, 'resil_lead_13'),
    ('resilience_30',  'scaffold_16', 563722, 'resil_bslmm_16'),
    ('resilience_30',  'scaffold_2',  71223,  'resil_bslmm_2'),
    ('AUC_20',         'scaffold_16', 421208, 'auc20_lead_16'),
    ('AUC_30',         'scaffold_4',  508257, 'auc30_lead_4'),
    ('cu_dose_slope',  'scaffold_3',  546068, 'cu_slope_lead_3'),
    ('IC50_est',       'scaffold_16', 122361, 'ic50_lead_16'),
]

WINDOW_KB = 150
PRIOR_SD = 0.2          # prior SD of causal NCP (scale-free, in z units)
W = PRIOR_SD ** 2


def z_abf(beta, se):
    """Wakefield approximate Bayes factor (log scale), scale-free in z = beta/se.

    Prior on the causal effect NCP: N(0, W).  ABF in z-space:
        logABF = 0.5*log(1/(1+W)) + 0.5*(W/(1+W))*z^2
    """
    if se is None or se <= 0 or beta is None:
        return -math.inf
    z = beta / se
    r = W / (1.0 + W)
    return 0.5 * math.log(1.0 - r) + 0.5 * r * z * z


def load_genes(path):
    with open(path) as f:
        d = json.load(f)
    by_scaf = {}
    for v in d.values():
        by_scaf.setdefault(v['scaf'], []).append(v)
    for k in by_scaf:
        by_scaf[k].sort(key=lambda x: x['start'])
    return d, by_scaf


def nearest_gene(by_scaf, scaf, pos):
    for g in by_scaf.get(scaf, []):
        if g['start'] <= pos <= g['end']:
            return g['gene'], 0, 'inside'
    hits = [g for g in by_scaf.get(scaf, []) if g['start'] > pos]
    if hits:
        g = min(hits, key=lambda x: x['start'] - pos)
        return g['gene'], g['start'] - pos, 'near'
    if by_scaf.get(scaf):
        g = max(by_scaf[scaf], key=lambda x: x['end'])
        return g['gene'], pos - g['end'], 'near'
    return None, None, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--assoc', default=ASSOC)
    ap.add_argument('--index', default=INDEX)
    ap.add_argument('--outdir', default=TIERE)
    ap.add_argument('--window-kb', type=int, default=WINDOW_KB)
    ap.add_argument('--cand-p', type=float, default=1e-3,
                    help='candidate variants must have p_wald < this (null markers otherwise '
                         'dilute the credible set to thousands of unassociated SNPs)')
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    genes, by_scaf = load_genes(args.index)
    window = args.window_kb * 1000

    rows = []
    cache = {}
    for trait, scaf, pos, label in ANCHORS:
        key = f'gwas_{trait}_assoc.csv.gz'
        if key not in cache:
            p = os.path.join(args.assoc, key)
            if not os.path.exists(p):
                print(f'  SKIP {trait}: no assoc file {p}', file=sys.stderr)
                continue
            df = pd.read_csv(p)
            df['scaf'] = 'scaffold_' + df['chr'].astype(str)
            cache[key] = df
        df = cache[key]
        win = df[(df['scaf'] == scaf) & (df['ps'] >= pos - window) & (df['ps'] <= pos + window)].copy()
        n_window = len(win)
        if n_window == 0:
            print(f'  SKIP {label}: no SNPs in window', file=sys.stderr)
            continue
        cand = win[win['p_wald'] < args.cand_p].copy()
        if cand.empty:
            print(f'  SKIP {label}: no candidates with p<{args.cand_p}', file=sys.stderr)
            continue
        win = cand.reset_index(drop=True)
        win['logABF'] = [z_abf(b, s) for b, s in zip(win['beta'], win['se'])]
        win = win.sort_values('logABF', ascending=False).reset_index(drop=True)
        lmax = float(win['logABF'].max())
        finite = win['logABF'][win['logABF'] != -math.inf]
        if finite.empty:
            print(f'  SKIP {label}: no finite ABF', file=sys.stderr)
            continue
        lse = float(lmax) if len(finite) == 1 else float(lmax) + math.log(sum(math.exp(x - lmax) for x in finite))
        win['pp'] = [math.exp(x - lse) if x != -math.inf else 0.0 for x in win['logABF']]
        win['cump'] = win['pp'].cumsum().round(8)

        gid, dist, rel = nearest_gene(by_scaf, scaf, pos)
        g = genes.get(gid, {}) if gid else {}
        lead = win.iloc[0]
        for level, thresh in (('90%', 0.90), ('95%', 0.95), ('99%', 0.99)):
            n_cs = max(1, int((win['cump'] <= thresh).sum()) + (0 if win.iloc[0]['cump'] <= thresh else 1))
            n_cs = min(n_cs, len(win))
            cs_rows = win.head(n_cs)
            lead_cs = cs_rows.iloc[0]
            low95, high95 = lead_cs['beta'] - 1.96 * lead_cs['se'], lead_cs['beta'] + 1.96 * lead_cs['se']
            rows.append({
                'label': label, 'trait': trait, 'scaffold': scaf, 'anchor_pos': pos,
                'n_snps_window': n_window, 'n_cand': len(win),
                'credset': level, 'n_snps_credset': n_cs,
                'lead_rs': lead['rs'], 'lead_pos': int(lead['ps']),
                'lead_pp': round(float(lead['pp']), 4), 'lead_af': float(lead['af']),
                'lead_beta': float(lead['beta']), 'lead_se': float(lead['se']),
                'beta_95ci_lo': round(float(low95), 4), 'beta_95ci_hi': round(float(high95), 4),
                'lead_p_wald': float(lead['p_wald']),
                'cs_max_beta': float(cs_rows['beta'].abs().max()),
                'cs_snp_min_pp': float(lead_cs['pp']),
                'nearest_gene': gid, 'gene_dist': dist, 'gene_rel': rel,
                'product': g.get('product'),
                'go': ';'.join(g.get('go', [])),
                'rare_driven': bool(lead['af'] < 0.02),
            })
        print(f'  {label}: {n_window} SNPs in {args.window_kb}kb, {len(win)} candidates, '
              f'99% CS n={[r["n_snps_credset"] for r in rows if r["label"]==label and r["credset"]=="99%"][0]}, '
              f'lead {win.iloc[0]["rs"]} pp={float(win.iloc[0]["pp"]):.3f} af={float(win.iloc[0]["af"]):.3f} beta={float(win.iloc[0]["beta"]):.2f}',
              file=sys.stderr)

    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(args.outdir, 'tierE_credible_sets.csv'), index=False)
    print(f'Wrote {len(rows)} rows -> {os.path.join(args.outdir, "tierE_credible_sets.csv")}', file=sys.stderr)


if __name__ == '__main__':
    main()
