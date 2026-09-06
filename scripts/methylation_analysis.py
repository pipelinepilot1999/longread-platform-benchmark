#!/usr/bin/env python3
"""Arm C analysis from modkit bedMethyl files (5mC, CpG, strands combined).

Outputs (all restricted to the analysed chromosome):
  methylation_global.tsv        histogram of per-site methylation + bimodality summary
  methylation_by_context.tsv    methylation by genomic context (islands, shores, promoters, gene bodies,
                                intergenic; promoter vs non-promoter islands; islands by distance to TSS)
  asm_sites.tsv                 per-CpG haplotype comparison (Fisher exact, BH-FDR); significant sites
  asm_sites_all.tsv.gz          every tested site
  asm_regions.tsv               merged runs of significant sites, annotated with known imprinted DMRs
  imprinting_control.tsv        the positive control: known imprinted DMRs unphased and per haplotype
  methylation_site_counts.tsv   how many sites survive each coverage threshold
  figures/fig07..fig10
"""
import argparse, gzip, os, sys
import numpy as np, pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
os.environ.setdefault('MPLCONFIGDIR', '/tmp/mpl')
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
for k in ['unphased', 'hp1', 'hp2', 'cpg_islands', 'gtf', 'imprinted', 'fai', 'chrom']:
    ap.add_argument('--' + k.replace('_', '-'), required=True)
ap.add_argument('--min-cov', type=int, default=10); ap.add_argument('--min-cov-hap', type=int, default=5)
ap.add_argument('--outdir', default='.')
a = ap.parse_args()
os.makedirs(os.path.join(a.outdir, 'figures'), exist_ok=True)
CH = a.chrom

def read_bedmethyl(path):
    cols = ['chrom', 'start', 'end', 'mod', 'score', 'strand', 's2', 'e2', 'color', 'depth', 'pct', 'nmod', 'ncanon', 'nother', 'ndel', 'nfail', 'ndiff', 'nnocall']
    df = pd.read_csv(path, sep='\t', header=None, names=cols, usecols=['chrom', 'start', 'mod', 'depth', 'pct', 'nmod', 'ncanon'])
    df = df[(df.chrom == CH) & (df['mod'] == 'm')].copy()
    df['frac'] = df.pct / 100.0
    return df[['start', 'depth', 'frac', 'nmod', 'ncanon']].set_index('start')

un = read_bedmethyl(a.unphased); h1 = read_bedmethyl(a.hp1); h2 = read_bedmethyl(a.hp2)
counts = []
for thr in [1, 3, 5, 10, 15, 20]:
    counts.append({'threshold': thr, 'unphased_sites': int((un.depth >= thr).sum()), 'hp1_sites': int((h1.depth >= thr).sum()),
                   'hp2_sites': int((h2.depth >= thr).sum()), 'sites_both_haplotypes': int(((h1.depth >= thr) & (h2.reindex(h1.index).depth >= thr)).sum())})
pd.DataFrame(counts).to_csv(os.path.join(a.outdir, 'methylation_site_counts.tsv'), sep='\t', index=False)

# ---------------- global distribution
u = un[un.depth >= a.min_cov]
hist, edges = np.histogram(u.frac, bins=np.linspace(0, 1, 21))
g = pd.DataFrame({'bin_start': edges[:-1], 'bin_end': edges[1:], 'n_sites': hist, 'fraction': hist / hist.sum()})
low = float((u.frac < 0.2).mean()); high = float((u.frac > 0.8).mean()); mid = float(((u.frac >= 0.4) & (u.frac <= 0.6)).mean())
g.attrs = {}
summary = pd.DataFrame([{'n_sites': len(u), 'min_cov': a.min_cov, 'mean_methylation': u.frac.mean(), 'median_methylation': u.frac.median(),
                         'fraction_lt_0.2': low, 'fraction_gt_0.8': high, 'fraction_0.4_0.6': mid,
                         'bimodal': (low + high) > 0.6 and mid < 0.15}])
with open(os.path.join(a.outdir, 'methylation_global.tsv'), 'w') as fh:
    summary.to_csv(fh, sep='\t', index=False); fh.write('\n'); g.to_csv(fh, sep='\t', index=False)
fig, ax = plt.subplots(figsize=(6, 4))
ax.hist(u.frac, bins=50, color='#4C72B0'); ax.set_xlabel('CpG methylation fraction (5mC, strands combined)'); ax.set_ylabel('CpG sites')
ax.set_title(f'{CH}: {len(u):,} CpGs with coverage >= {a.min_cov}\n<0.2: {low:.1%}   >0.8: {high:.1%}   0.4-0.6: {mid:.1%}')
fig.tight_layout(); fig.savefig(os.path.join(a.outdir, 'figures/fig07_global_methylation_distribution.png'), dpi=150); plt.close(fig)

# ---------------- annotation
isl = pd.read_csv(a.cpg_islands, sep='\t', header=None, usecols=[1, 2, 3], names=['chrom', 'start', 'end'])
isl = isl[isl.chrom == CH].sort_values('start').reset_index(drop=True)
genes = []
with gzip.open(a.gtf, 'rt') as fh:
    for line in fh:
        if line.startswith('#'): continue
        f = line.rstrip('\n').split('\t')
        if f[0] != CH or f[2] != 'gene' or 'gene_type "protein_coding"' not in f[8]: continue
        name = f[8].split('gene_name "')[1].split('"')[0]
        s, e, strand = int(f[3]) - 1, int(f[4]), f[6]
        genes.append({'name': name, 'start': s, 'end': e, 'strand': strand, 'tss': s if strand == '+' else e - 1})
genes = pd.DataFrame(genes)

def interval_mask(pos, starts, ends):
    """boolean: does each pos fall inside any [start,end) interval (intervals sorted, may overlap)."""
    starts = np.asarray(starts); ends = np.asarray(ends)
    order = np.argsort(starts); starts, ends = starts[order], ends[order]
    cummax_end = np.maximum.accumulate(ends)
    i = np.searchsorted(starts, pos, side='right') - 1
    ok = i >= 0
    res = np.zeros(len(pos), bool)
    res[ok] = pos[ok] < cummax_end[i[ok]]
    # cummax handles overlaps conservatively (a pos inside an earlier, longer interval is still caught)
    return res

pos = u.index.values
in_island = interval_mask(pos, isl.start, isl.end)
in_shore = interval_mask(pos, np.maximum(isl.start - 2000, 0), isl.end + 2000) & ~in_island
prom = pd.DataFrame({'start': genes.tss - 1000, 'end': genes.tss + 1000})
in_prom = interval_mask(pos, prom.start, prom.end)
in_body = interval_mask(pos, genes.start, genes.end)
ctx = np.where(in_island, 'cpg_island', np.where(in_shore, 'cpg_shore', np.where(in_prom, 'promoter_non_island', np.where(in_body, 'gene_body', 'intergenic'))))
ctx_df = pd.DataFrame({'pos': pos, 'frac': u.frac.values, 'context': ctx})
rows = []
for c, sub in ctx_df.groupby('context'):
    rows.append({'section': 'context_exclusive', 'category': c, 'n_sites': len(sub), 'mean_methylation': sub.frac.mean(), 'median_methylation': sub.frac.median(),
                 'fraction_lt_0.2': (sub.frac < 0.2).mean(), 'fraction_gt_0.8': (sub.frac > 0.8).mean()})
# promoter islands vs other islands (site-level)
isl_prom = interval_mask(pos, isl.start, isl.end) & in_prom
for name, m in [('promoter_cpg_island', in_island & in_prom), ('non_promoter_cpg_island', in_island & ~in_prom), ('promoter_all', in_prom), ('non_promoter_all', ~in_prom)]:
    sub = u.frac.values[m]
    rows.append({'section': 'promoter_islands', 'category': name, 'n_sites': int(m.sum()), 'mean_methylation': sub.mean(), 'median_methylation': np.median(sub),
                 'fraction_lt_0.2': (sub < 0.2).mean(), 'fraction_gt_0.8': (sub > 0.8).mean()})
# islands by distance to nearest TSS (island-level mean methylation)
tss = np.sort(genes.tss.values)
isl_mean = []
for s, e in zip(isl.start, isl.end):
    sel = u.frac.values[(pos >= s) & (pos < e)]
    if len(sel) < 5: continue
    mid = (s + e) // 2; j = np.searchsorted(tss, mid); d = min(abs(tss[max(j - 1, 0)] - mid), abs(tss[min(j, len(tss) - 1)] - mid))
    isl_mean.append({'start': s, 'end': e, 'n_cpg': len(sel), 'mean': sel.mean(), 'dist_tss': d})
isl_mean = pd.DataFrame(isl_mean)
bins = [(-1, 1000, '<=1kb'), (1000, 5000, '1-5kb'), (5000, 20000, '5-20kb'), (20000, 100000, '20-100kb'), (100000, 10**9, '>100kb')]
for lo, hi, nm in bins:
    sub = isl_mean[(isl_mean.dist_tss > lo) & (isl_mean.dist_tss <= hi)]
    rows.append({'section': 'islands_by_tss_distance', 'category': nm, 'n_sites': len(sub), 'mean_methylation': sub['mean'].mean() if len(sub) else np.nan,
                 'median_methylation': sub['mean'].median() if len(sub) else np.nan, 'fraction_lt_0.2': (sub['mean'] < 0.2).mean() if len(sub) else np.nan,
                 'fraction_gt_0.8': (sub['mean'] > 0.8).mean() if len(sub) else np.nan})
pd.DataFrame(rows).to_csv(os.path.join(a.outdir, 'methylation_by_context.tsv'), sep='\t', index=False)
order = ['cpg_island', 'cpg_shore', 'promoter_non_island', 'gene_body', 'intergenic']
fig, axes = plt.subplots(1, 2, figsize=(11, 4), gridspec_kw={'width_ratios': [3, 2]})
axes[0].violinplot([ctx_df.frac[ctx_df.context == c].values for c in order], showmedians=True)
axes[0].set_xticks(range(1, 6)); axes[0].set_xticklabels([f"{c}\n(n={int((ctx_df.context==c).sum()):,})" for c in order], fontsize=8)
axes[0].set_ylabel('CpG methylation fraction'); axes[0].set_title('Methylation by genomic context (exclusive assignment)')
labels = [b[2] for b in bins]; vals = [isl_mean[(isl_mean.dist_tss > lo) & (isl_mean.dist_tss <= hi)]['mean'].values for lo, hi, _ in bins]
axes[1].boxplot(vals, labels=labels); axes[1].set_xlabel('CpG island distance to nearest protein-coding TSS'); axes[1].set_ylabel('island mean methylation')
axes[1].set_title('CpG islands by distance to TSS'); axes[1].tick_params(axis='x', labelsize=8)
fig.tight_layout(); fig.savefig(os.path.join(a.outdir, 'figures/fig08_methylation_by_context.png'), dpi=150); plt.close(fig)

# ---------------- allele-specific methylation
j = h1.join(h2, lsuffix='_1', rsuffix='_2', how='inner')
j = j[(j.depth_1 >= a.min_cov_hap) & (j.depth_2 >= a.min_cov_hap)].copy()
def fisher(r):
    return fisher_exact([[r.nmod_1, r.ncanon_1], [r.nmod_2, r.ncanon_2]])[1]
# vectorised-ish: only test sites where a difference is possible; identical rows get p=1
j['delta'] = j.frac_1 - j.frac_2
j['p'] = 1.0
mask = (j.nmod_1 + j.ncanon_1 > 0) & (j.nmod_2 + j.ncanon_2 > 0)
j.loc[mask, 'p'] = [fisher_exact([[int(m1), int(c1)], [int(m2), int(c2)]])[1] for m1, c1, m2, c2 in zip(j.nmod_1[mask], j.ncanon_1[mask], j.nmod_2[mask], j.ncanon_2[mask])]
j['q'] = multipletests(j.p, method='fdr_bh')[1] if len(j) else []
j['significant'] = (j.q < 0.05) & (j.delta.abs() >= 0.5)
j['pos'] = j.index
imp = pd.read_csv(a.imprinted, sep='\t', comment='#', header=None, names=['chrom', 'start', 'end', 'name', 'origin', 'source', 'class'])
imp = imp[imp.chrom == CH]
known = imp[imp['class'].isin(['primary_control', 'known_other'])]
def imprinted_label(p):
    hit = known[(known.start <= p) & (p < known.end)]
    return ','.join(hit.name) if len(hit) else ''
j['known_imprinted_dmr'] = [imprinted_label(p) for p in j.pos]
cols = ['pos', 'depth_1', 'frac_1', 'depth_2', 'frac_2', 'delta', 'p', 'q', 'significant', 'known_imprinted_dmr']
j[cols].to_csv(os.path.join(a.outdir, 'asm_sites_all.tsv.gz'), sep='\t', index=False, float_format='%.4g')
j[j.significant][cols].to_csv(os.path.join(a.outdir, 'asm_sites.tsv'), sep='\t', index=False, float_format='%.4g')

# merge significant sites into regions (same sign, gap <= 500 bp, >= 5 CpGs)
sig = j[j.significant].sort_values('pos')
regions = []
cur = None
for p, d in zip(sig.pos, sig.delta):
    sgn = 1 if d > 0 else -1
    if cur and sgn == cur['sign'] and p - cur['end'] <= 500:
        cur['end'] = p; cur['n'] += 1; cur['deltas'].append(d)
    else:
        if cur: regions.append(cur)
        cur = {'start': p, 'end': p, 'n': 1, 'sign': sgn, 'deltas': [d]}
if cur: regions.append(cur)
reg_rows = []
for r in regions:
    if r['n'] < 5: continue
    hit = known[(known.start < r['end'] + 1) & (r['start'] < known.end)]
    span = (r['start'], r['end'] + 1)
    near = genes[(genes.start < span[1] + 5000) & (genes.end > span[0] - 5000)]
    tested = j[(j.pos >= span[0]) & (j.pos < span[1])]
    reg_rows.append({'chrom': CH, 'start': span[0], 'end': span[1], 'n_significant_cpg': r['n'], 'n_tested_cpg': len(tested),
                     'mean_delta_hp1_minus_hp2': np.mean(r['deltas']), 'mean_frac_hp1': tested.frac_1.mean(), 'mean_frac_hp2': tested.frac_2.mean(),
                     'overlaps_known_imprinted_dmr': ','.join(hit.name) if len(hit) else '', 'genes_within_5kb': ','.join(near.name.head(5)),
                     'classification': 'known_imprinted_control' if len(hit) else 'candidate_ASM'})
reg = pd.DataFrame(reg_rows).sort_values('n_significant_cpg', ascending=False) if reg_rows else pd.DataFrame(columns=['chrom'])
reg.to_csv(os.path.join(a.outdir, 'asm_regions.tsv'), sep='\t', index=False, float_format='%.4f')

# ---------------- imprinting control
ctrl = []
for _, d in imp.iterrows():
    su = un[(un.index >= d.start) & (un.index < d.end) & (un.depth >= a.min_cov)]
    sj = j[(j.pos >= d.start) & (j.pos < d.end)]
    w = lambda s, f, c: float(np.average(s[f], weights=s[c])) if len(s) else np.nan
    m_un = w(su, 'frac', 'depth'); m1 = w(sj, 'frac_1', 'depth_1'); m2 = w(sj, 'frac_2', 'depth_2')
    sep = abs(m1 - m2) if len(sj) else np.nan
    ok = (len(su) >= 5 and 0.25 <= m_un <= 0.75) and (len(sj) >= 5 and sep >= 0.5 and max(m1, m2) >= 0.7 and min(m1, m2) <= 0.3)
    ctrl.append({'name': d['name'], 'class': d['class'], 'origin_of_methylated_allele': d.origin, 'start': d.start, 'end': d.end,
                 'n_cpg_unphased': len(su), 'unphased_methylation': m_un, 'n_cpg_both_haplotypes': len(sj),
                 'hp1_methylation': m1, 'hp2_methylation': m2, 'haplotype_separation': sep,
                 'methylated_haplotype': ('hp1' if m1 > m2 else 'hp2') if len(sj) else 'NA',
                 'fraction_sites_significant': sj.significant.mean() if len(sj) else np.nan,
                 'control_pass': bool(ok)})
ctrl = pd.DataFrame(ctrl)
ctrl.to_csv(os.path.join(a.outdir, 'imprinting_control.tsv'), sep='\t', index=False, float_format='%.4f')
prim = ctrl[ctrl['class'] == 'primary_control']
print('IMPRINTING CONTROL:', int(prim.control_pass.sum()), 'of', len(prim), 'primary-control DMRs pass')
print(prim[['name', 'n_cpg_unphased', 'unphased_methylation', 'hp1_methylation', 'hp2_methylation', 'haplotype_separation', 'control_pass']].to_string())

# fig 9: GNAS locus
lo, hi = 58_800_000, 58_920_000
fig, ax = plt.subplots(figsize=(12, 4.5))
su = un[(un.index >= lo) & (un.index < hi) & (un.depth >= a.min_cov)]
sj = j[(j.pos >= lo) & (j.pos < hi)]
ax.scatter(su.index, su.frac, s=4, c='grey', alpha=0.5, label='unphased')
ax.scatter(sj.pos, sj.frac_1, s=6, c='#D55E00', alpha=0.7, label='haplotype 1')
ax.scatter(sj.pos, sj.frac_2, s=6, c='#0072B2', alpha=0.7, label='haplotype 2')
for _, d in imp[(imp['class'] == 'primary_control') & (imp.start > lo) & (imp.end < hi)].iterrows():
    ax.axvspan(d.start, d.end, color='gold', alpha=0.25); ax.text((d.start + d.end) / 2, 1.04, d['name'].split(',')[0], ha='center', fontsize=7, rotation=0)
for _, gname in genes[(genes.start < hi) & (genes.end > lo)].iterrows():
    ax.plot([gname.start, gname.end], [-0.06, -0.06], lw=3, c='black'); ax.text(gname.start, -0.12, gname['name'], fontsize=7)
ax.set_ylim(-0.15, 1.1); ax.set_xlim(lo, hi); ax.set_xlabel(f'{CH} position (bp)'); ax.set_ylabel('CpG methylation')
ax.set_title('GNAS imprinted locus: per-CpG methylation, unphased vs haplotype-partitioned (known DMRs shaded)'); ax.legend(loc='center left', fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(a.outdir, 'figures/fig09_gnas_imprinted_locus.png'), dpi=150); plt.close(fig)

# fig 10: volcano
fig, ax = plt.subplots(figsize=(6.5, 5))
ql = -np.log10(np.clip(j.q.values, 1e-300, 1))
ax.scatter(j.delta, ql, s=2, c='lightgrey', label='tested CpG')
m = j.known_imprinted_dmr != ''
ax.scatter(j.delta[m], ql[m], s=6, c='#D55E00', label='in known imprinted DMR')
ax.axhline(-np.log10(0.05), ls='--', c='k', lw=0.8); ax.axvline(0.5, ls='--', c='k', lw=0.8); ax.axvline(-0.5, ls='--', c='k', lw=0.8)
ax.set_xlabel('methylation difference (haplotype 1 - haplotype 2)'); ax.set_ylabel('-log10 BH-adjusted p (Fisher exact)')
ax.set_title(f'Allele-specific methylation across {CH}: {int(j.significant.sum()):,} significant of {len(j):,} tested CpGs'); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(a.outdir, 'figures/fig10_asm_volcano.png'), dpi=150); plt.close(fig)
print('ASM sites tested', len(j), 'significant', int(j.significant.sum()), 'regions', len(reg))
