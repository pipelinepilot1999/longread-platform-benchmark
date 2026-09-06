#!/usr/bin/env python3
"""Cross-platform and cross-arm comparison (run after both platform runs).

Reads results/ont and results/hifi (plus the Truvari platform-concordance directory) and writes:
  qc_summary.tsv, sv_benchmark.tsv, sv_concordance.tsv, smallvar_benchmark.tsv, phasing_stats.tsv  (both platforms merged)
  sv_platform_by_size.tsv      recall by SV type and size bin, ONT vs HiFi, with the crossover
  platform_concordance.tsv     calls found by both platforms / ONT-only / HiFi-only, cross-referenced to truth
  cross_arm.tsv                section 9 tests with effect sizes
  figures/fig01..fig06
"""
import argparse, glob, gzip, os, sys, json
import numpy as np, pandas as pd, pysam
from scipy import stats
os.environ.setdefault('MPLCONFIGDIR', '/tmp/mpl')
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(__file__))
from sv_common import read_vcf, SIZE_BINS

ap = argparse.ArgumentParser()
ap.add_argument('--ont', required=True); ap.add_argument('--hifi', required=True); ap.add_argument('--concordance', required=True)
ap.add_argument('--chrom', default='chr20'); ap.add_argument('--chrom-len', type=int, default=64444167)
ap.add_argument('--best', default='sniffles,sniffles', help='best caller ont,hifi'); ap.add_argument('--outdir', default='.')
ap.add_argument('--seed', type=int, default=11)
a = ap.parse_args()
os.makedirs(os.path.join(a.outdir, 'figures'), exist_ok=True)
rng = np.random.default_rng(a.seed)
P = {'ont': a.ont, 'hifi': a.hifi}
COL = {'ont': '#D55E00', 'hifi': '#0072B2'}
best = dict(zip(['ont', 'hifi'], a.best.split(',')))
def rd(p, name): return pd.read_csv(os.path.join(P[p], name), sep='\t')

# ---------------- merged tables
for name in ['qc_summary', 'sv_benchmark', 'sv_concordance', 'smallvar_benchmark', 'phasing_stats']:
    pd.concat([rd(p, f'{name}.{p}.tsv') for p in P]).to_csv(os.path.join(a.outdir, f'{name}.tsv'), sep='\t', index=False)
qc = pd.read_csv(os.path.join(a.outdir, 'qc_summary.tsv'), sep='\t').set_index('platform')
svb = pd.read_csv(os.path.join(a.outdir, 'sv_benchmark.tsv'), sep='\t')
smv = pd.read_csv(os.path.join(a.outdir, 'smallvar_benchmark.tsv'), sep='\t')
ph = pd.read_csv(os.path.join(a.outdir, 'phasing_stats.tsv'), sep='\t').set_index('platform')

# ---------------- fig 1: read lengths + depth profile
fig, axes = plt.subplots(1, 2, figsize=(12, 4))
for p in P:
    rl = np.array([int(x) for x in gzip.open(os.path.join(P[p], 'qc', f'{p}.read_lengths.txt.gz'), 'rt')])
    axes[0].hist(rl, bins=np.logspace(2, 6, 80), histtype='step', lw=1.5, color=COL[p], label=f"{p.upper()} (N50 {int(qc.loc[p,'read_n50']):,} bp)", weights=rl)
    reg = pd.read_csv(os.path.join(P[p], 'qc', f'{p}.regions.bed.gz'), sep='\t', header=None, names=['c', 's', 'e', 'd'])
    reg = reg[reg.c == a.chrom]
    axes[1].plot(reg.s / 1e6, reg.d, lw=0.6, color=COL[p], alpha=0.8, label=f"{p.upper()} (mean {qc.loc[p,'mean_depth_target_chrom']:.1f}x)")
axes[0].set_xscale('log'); axes[0].set_xlabel('read length (bp)'); axes[0].set_ylabel('bases (weighted histogram)'); axes[0].legend(); axes[0].set_title('Read length distribution (analysed reads)')
axes[1].set_ylim(0, max(float(qc.loc[p, 'mean_depth_target_chrom']) for p in P) * 2.5); axes[1].set_xlabel(f'{a.chrom} position (Mb)'); axes[1].set_ylabel('depth (10 kb bins)'); axes[1].legend(); axes[1].set_title('Depth profile after coverage matching')
fig.tight_layout(); fig.savefig(os.path.join(a.outdir, 'figures/fig01_read_length_and_depth.png'), dpi=150); plt.close(fig)

# ---------------- fig 2: SV precision/recall by caller, platform, type, regime
regimes = ['permissive', 'default', 'strict']
fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), sharey=True)
for ax, t in zip(axes, ['ALL', 'DEL', 'INS']):
    sub = svb[(svb.svtype == t) & (svb.size_bin == 'ALL')]
    for p in P:
        for caller, mk in [('sniffles', 'o'), ('cutesv', 's')]:
            s = sub[(sub.platform == p) & (sub.caller == caller)].set_index('regime').loc[regimes]
            ax.plot(s.recall, s.precision, '-', color=COL[p], alpha=0.4, lw=1)
            for rg, ms in zip(regimes, [4, 7, 10]):
                ax.plot(s.loc[rg, 'recall'], s.loc[rg, 'precision'], mk, color=COL[p], ms=ms, mfc='none' if caller == 'cutesv' else COL[p], label=f'{p.upper()} {caller}' if rg == 'default' else None)
    ax.set_title(f'SV type: {t}'); ax.set_xlabel('recall'); ax.grid(alpha=0.3)
axes[0].set_ylabel('precision'); axes[0].legend(fontsize=7, loc='lower left')
fig.suptitle('Truvari precision/recall vs GIAB v5.0q (benchmark BED). Marker size: permissive < default < strict regime; filled = Sniffles2, open = cuteSV')
fig.tight_layout(); fig.savefig(os.path.join(a.outdir, 'figures/fig02_sv_precision_recall.png'), dpi=150); plt.close(fig)

# ---------------- fig 4 + sv_platform_by_size: recall by size bin ONT vs HiFi (best caller, default regime)
rows = []
for t in ['DEL', 'INS', 'ALL']:
    for b in [x[2] for x in SIZE_BINS] + ['ALL']:
        r = {'svtype': t, 'size_bin': b}
        for p in P:
            s = svb[(svb.platform == p) & (svb.caller == best[p]) & (svb.regime == 'default') & (svb.svtype == t) & (svb.size_bin == b)]
            if len(s):
                s = s.iloc[0]
                r.update({f'{p}_caller': best[p], f'{p}_tp': int(s.tp_base), f'{p}_fn': int(s.fn), f'{p}_fp': int(s.fp), f'{p}_recall': float(s.recall), f'{p}_precision': float(s.precision), f'{p}_f1': float(s.f1)})
        r['recall_diff_ont_minus_hifi'] = r.get('ont_recall', np.nan) - r.get('hifi_recall', np.nan)
        rows.append(r)
sz = pd.DataFrame(rows)
def crossover(t):
    s = sz[(sz.svtype == t) & (sz.size_bin != 'ALL')].reset_index(drop=True)
    d = s.recall_diff_ont_minus_hifi.values; bins = s.size_bin.tolist()
    for i in range(1, len(d)):
        if np.sign(d[i]) != np.sign(d[i - 1]) and not (np.isnan(d[i]) or np.isnan(d[i - 1])) and d[i] != 0:
            return f'between {bins[i-1]} and {bins[i]}', i
    return ('no crossover (ONT >= HiFi in every bin)' if np.nanmin(d) >= 0 else 'no crossover (HiFi >= ONT in every bin)' if np.nanmax(d) <= 0 else 'no clean crossover'), None
sz['crossover'] = [crossover(t)[0] if b == 'ALL' else '' for t, b in zip(sz.svtype, sz.size_bin)]
sz.to_csv(os.path.join(a.outdir, 'sv_platform_by_size.tsv'), sep='\t', index=False, float_format='%.4f')
fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
labels = [x[2] for x in SIZE_BINS]
for ax, t in zip(axes, ['DEL', 'INS']):
    s = sz[(sz.svtype == t) & (sz.size_bin != 'ALL')]
    x = np.arange(len(labels)); w = 0.38
    for k, p in enumerate(P):
        ax.bar(x + (k - 0.5) * w, s[f'{p}_recall'], w, color=COL[p], label=f'{p.upper()} ({best[p]})')
        for xi, (tp, fn) in enumerate(zip(s[f'{p}_tp'], s[f'{p}_fn'])):
            ax.text(xi + (k - 0.5) * w, 0.02, f'{tp}/{tp+fn}', ha='center', fontsize=6, rotation=90, color='white')
    cx, ci = crossover(t)
    if ci is not None:
        ax.axvline(ci - 0.5, ls='--', color='k'); ax.text(ci - 0.5, 1.02, 'crossover', ha='center', fontsize=8)
    ax.set_xticks(x); ax.set_xticklabels(labels); ax.set_title(f'{t}: recall by size bin ({cx})', fontsize=9); ax.set_xlabel('SV size'); ax.set_ylim(0, 1.08)
axes[0].set_ylabel('recall (Truvari default regime, benchmark BED)'); axes[0].legend(fontsize=8)
fig.suptitle('SV recall by size: ONT vs PacBio HiFi at matched coverage (TP/total truth in bars)')
fig.tight_layout(); fig.savefig(os.path.join(a.outdir, 'figures/fig04_sv_recall_by_size_crossover.png'), dpi=150); plt.close(fig)

# ---------------- platform concordance (Truvari ONT-best vs HiFi-best inside the BED) + fig 3
calls = {p: pd.read_csv(os.path.join(P[p], f'sv_calls_annotated.{p}.tsv.gz'), sep='\t') for p in P}
def truth_of(p, caller, rid):
    c = calls[p]; s = c[(c.caller == caller) & (c.id == rid)]
    return s.status_default.iloc[0] if len(s) and isinstance(s.status_default.iloc[0], str) else 'not_benchmarked'
cd = a.concordance
groups = [('both_platforms', 'ont', read_vcf(os.path.join(cd, 'tp-base.vcf.gz'))), ('both_platforms', 'hifi', read_vcf(os.path.join(cd, 'tp-comp.vcf.gz'))),
          ('ont_only', 'ont', read_vcf(os.path.join(cd, 'fn.vcf.gz'))), ('hifi_only', 'hifi', read_vcf(os.path.join(cd, 'fp.vcf.gz')))]
pc_rows = []
for cat, p, recs in groups:
    for t in ['ALL', 'DEL', 'INS']:
        sel = [r for r in recs if t == 'ALL' or r['svtype'] == t]
        st = pd.Series([truth_of(p, best[p], r['id']) for r in sel]).value_counts()
        tp, fp = int(st.get('TP', 0)), int(st.get('FP', 0))
        pc_rows.append({'category': cat, 'representative_platform': p, 'caller': best[p], 'svtype': t, 'n_calls': len(sel), 'truth_TP': tp, 'truth_FP': fp,
                        'not_benchmarked': int(st.get('not_benchmarked', 0)), 'fp_fraction_of_benchmarked': tp + fp and fp / (tp + fp)})
pc = pd.DataFrame(pc_rows); pc.to_csv(os.path.join(a.outdir, 'platform_concordance.tsv'), sep='\t', index=False, float_format='%.4f')
cc = pd.read_csv(os.path.join(a.outdir, 'sv_concordance.tsv'), sep='\t')
fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
def stacked(ax, df, cats, title):
    df = df[df.svtype == 'ALL']
    xs = np.arange(len(cats))
    tp = [int(df[df.category == c].truth_TP.iloc[0]) for c in cats]; fp = [int(df[df.category == c].truth_FP.iloc[0]) for c in cats]
    nb = [int(df[df.category == c].not_benchmarked.iloc[0]) for c in cats]
    ax.bar(xs, tp, color='#2ca02c', label='TP vs truth'); ax.bar(xs, fp, bottom=tp, color='#d62728', label='FP vs truth'); ax.bar(xs, nb, bottom=np.add(tp, fp), color='lightgrey', label='outside benchmark')
    for x, (t, f) in enumerate(zip(tp, fp)):
        ax.text(x, t + f + 5, f'FP {f/(t+f):.0%}' if t + f else '', ha='center', fontsize=8)
    ax.set_xticks(xs); ax.set_xticklabels(cats, fontsize=8); ax.set_title(title, fontsize=9)
for ax, p in zip(axes[:2], P):
    d = cc[(cc.platform == p) & ((cc.category != 'shared') | (cc.representative_caller == 'sniffles'))]
    stacked(ax, d, ['shared', 'sniffles_only', 'cutesv_only'], f'{p.upper()}: Sniffles2 vs cuteSV (default regime, benchmark BED)')
d = pc[(pc.category != 'both_platforms') | (pc.representative_platform == 'ont')]
stacked(axes[2], d, ['both_platforms', 'ont_only', 'hifi_only'], f'Platforms: ONT {best["ont"]} vs HiFi {best["hifi"]}')
axes[0].set_ylabel('SV calls'); axes[0].legend(fontsize=7)
fig.tight_layout(); fig.savefig(os.path.join(a.outdir, 'figures/fig03_sv_concordance.png'), dpi=150); plt.close(fig)

# ---------------- fig 5: small variants with homopolymer stratification
strata = ['*', 'not_in_homopolymers_ge7', 'homopolymer_4to6', 'homopolymer_7to11', 'homopolymer_ge12', 'all_TR_and_homopolymers']
fig, axes = plt.subplots(1, 2, figsize=(13, 4.2), sharey=True)
for ax, t in zip(axes, ['SNP', 'INDEL']):
    x = np.arange(len(strata)); w = 0.2
    for k, (p, metric) in enumerate([(p, m) for p in P for m in ['recall', 'precision']]):
        vals = [smv[(smv.platform == p) & (smv.type == t) & (smv.stratum == s)][metric].iloc[0] if len(smv[(smv.platform == p) & (smv.type == t) & (smv.stratum == s)]) else np.nan for s in strata]
        ax.bar(x + (k - 1.5) * w, vals, w, color=COL[p], alpha=1 if metric == 'recall' else 0.45, label=f'{p.upper()} {metric}')
    ax.set_xticks(x); ax.set_xticklabels(['all', 'not in homopol≥ 7', 'homopol 4-6', 'homopol 7-11', 'homopol ≥ 12', 'any TR/homopol'], fontsize=8, rotation=15)
    ax.set_title(f'{t}: Clair3 vs GIAB v5.0q (hap.py, benchmark BED)'); ax.set_ylim(0, 1.05); ax.grid(axis='y', alpha=0.3)
axes[0].set_ylabel('metric'); axes[0].legend(fontsize=7, loc='lower left')
fig.tight_layout(); fig.savefig(os.path.join(a.outdir, 'figures/fig05_small_variants_homopolymer.png'), dpi=150); plt.close(fig)

# ---------------- fig 6: phase block length distributions
fig, ax = plt.subplots(figsize=(7, 4.2))
for p in P:
    bl = pd.read_csv(os.path.join(P[p], 'phasing', f'{p}.blocks.tsv'), sep='\t'); L = (bl['to'] - bl['from'] + 1).values
    ax.hist(L, bins=np.logspace(2, 8, 60), histtype='step', lw=1.6, color=COL[p], weights=L, label=f"{p.upper()}: {len(L):,} blocks, N50 {int(ph.loc[p,'block_n50_bp']):,} bp")
    ax.axvline(ph.loc[p, 'block_n50_bp'], color=COL[p], ls='--', lw=1)
ax.set_xscale('log'); ax.set_xlabel('phase block length (bp)'); ax.set_ylabel('bp in blocks (weighted)'); ax.legend(fontsize=8); ax.set_title('WhatsHap phase blocks on chr20 (dashed = N50)')
fig.tight_layout(); fig.savefig(os.path.join(a.outdir, 'figures/fig06_phase_block_lengths.png'), dpi=150); plt.close(fig)

# ---------------- cross-arm (section 9)
xa = []
def blocks_of(p):
    bl = pd.read_csv(os.path.join(P[p], 'phasing', f'{p}.blocks.tsv'), sep='\t'); return bl.sort_values('from')
def tp_svs(p):
    c = calls[p]; c = c[(c.caller == best[p]) & (c.status_default == 'TP')]; return c
for p in P:
    bl = blocks_of(p); sv = tp_svs(p)
    gaps = [(int(e), int(s)) for e, s in zip(bl['to'].values[:-1], bl['from'].values[1:]) if s > e]
    svpos = np.sort(sv.pos.values)
    def frac_gaps_with_sv(gaps, pad=1000):
        return np.mean([np.searchsorted(svpos, hi + pad) > np.searchsorted(svpos, lo - pad) for lo, hi in gaps]) if gaps else np.nan
    obs = frac_gaps_with_sv(gaps)
    null = []
    for _ in range(500):
        shift = rng.integers(0, a.chrom_len)
        null.append(frac_gaps_with_sv([((lo + shift) % a.chrom_len, (hi + shift) % a.chrom_len) for lo, hi in gaps if (hi + shift) % a.chrom_len > (lo + shift) % a.chrom_len]))
    null = np.array(null)
    xa.append({'platform': p, 'test': 'phase_block_gaps_containing_TP_SV_within_1kb', 'n_gaps': len(gaps), 'n_tp_sv': len(svpos), 'observed_fraction': obs,
               'null_mean_fraction': null.mean(), 'null_sd': null.std(), 'enrichment': obs / null.mean() if null.mean() else np.nan,
               'empirical_p': float((null >= obs).mean()), 'note': 'null = circularly shifted gap positions (500 shifts)'})

# haplotag quality vs heterozygosity (500 kb windows) using the haplotagged BAM and phased VCF
for p in P:
    bam = glob.glob(os.path.join(P[p], 'alignment', f'{p}.{a.chrom}.haplotagged.bam'))
    vcf = glob.glob(os.path.join(P[p], 'phasing', f'{p}.phased.vcf.gz'))
    if not bam or not vcf: continue
    W = 500_000; nw = a.chrom_len // W + 1
    tot = np.zeros(nw); un = np.zeros(nw); het = np.zeros(nw)
    with pysam.AlignmentFile(bam[0]) as bf:
        for r in bf.fetch(a.chrom):
            if r.is_secondary or r.is_supplementary or r.is_unmapped: continue
            w = r.reference_start // W; tot[w] += 1
            if not r.has_tag('HP'): un[w] += 1
    with pysam.VariantFile(vcf[0]) as vf:
        for rec in vf.fetch(a.chrom):
            gt = rec.samples[0]['GT']
            if gt and None not in gt and len(set(gt)) > 1: het[rec.pos // W] += 1
    ok = tot > 200
    rho, pv = stats.spearmanr(het[ok], un[ok] / tot[ok])
    xa.append({'platform': p, 'test': 'unassigned_read_fraction_vs_het_density_500kb', 'n_gaps': int(ok.sum()), 'n_tp_sv': '', 'observed_fraction': float(np.median(un[ok] / tot[ok])),
               'null_mean_fraction': '', 'null_sd': '', 'enrichment': rho, 'empirical_p': pv, 'note': 'enrichment column = Spearman rho between het variants per window and unassigned fraction; observed = median unassigned fraction'})
    pd.DataFrame({'window_start': np.arange(nw) * W, 'reads': tot, 'unassigned': un, 'het_variants': het}).to_csv(os.path.join(a.outdir, f'haplotag_windows.{p}.tsv'), sep='\t', index=False)

# methylation near SV breakpoints (ONT): windows +/-1 kb around TP SV breakpoints vs CpG-density-matched random windows
mfile = glob.glob(os.path.join(P['ont'], 'methylation', 'ont.cpg.unphased.bed.gz'))
if mfile:
    m = pd.read_csv(mfile[0], sep='\t', header=None, usecols=[0, 1, 3, 9, 10], names=['c', 's', 'mod', 'cov', 'pct'])
    m = m[(m.c == a.chrom) & (m['mod'] == 'm') & (m['cov'] >= 10)]; mp = m.s.values; mf = m.pct.values / 100
    sv = tp_svs('ont'); bp = np.unique(np.concatenate([sv.pos.values, sv.pos.values + np.where(sv.svtype == 'DEL', sv.svlen.values, 0)]))
    def win_stats(centers, pad=1000):
        out = []
        for c in centers:
            i, j = np.searchsorted(mp, c - pad), np.searchsorted(mp, c + pad)
            if j - i >= 5: out.append((j - i, mf[i:j].mean()))
        return np.array(out)
    obs = win_stats(bp)
    ctrl_c = rng.integers(1000, a.chrom_len - 1000, size=len(bp) * 20); ctrl = win_stats(ctrl_c)
    # match on CpG count decile
    if len(obs) and len(ctrl):
        q = np.quantile(obs[:, 0], np.linspace(0, 1, 6)); matched = []
        for k in range(5):
            o = obs[(obs[:, 0] >= q[k]) & (obs[:, 0] <= q[k + 1])]; c = ctrl[(ctrl[:, 0] >= q[k]) & (ctrl[:, 0] <= q[k + 1])]
            if len(o) and len(c): matched.append(c[rng.choice(len(c), size=min(len(c), len(o) * 5), replace=False)])
        matched = np.concatenate(matched)
        u, pv = stats.mannwhitneyu(obs[:, 1], matched[:, 1])
        xa.append({'platform': 'ont', 'test': 'methylation_within_1kb_of_TP_SV_breakpoints_vs_CpG_matched_random_windows', 'n_gaps': len(obs), 'n_tp_sv': len(matched),
                   'observed_fraction': obs[:, 1].mean(), 'null_mean_fraction': matched[:, 1].mean(), 'null_sd': matched[:, 1].std(),
                   'enrichment': obs[:, 1].mean() - matched[:, 1].mean(), 'empirical_p': pv, 'note': 'enrichment column = mean methylation difference (SV windows - control); Mann-Whitney p'})

# platform disagreement overlap: 100 kb windows with platform-unique SV calls vs windows with switch errors / block boundaries in one platform only
W = 100_000; nw = a.chrom_len // W + 1
sv_div = np.zeros(nw, bool)
for cat, p, recs in groups:
    if cat in ('ont_only', 'hifi_only'):
        for r in recs: sv_div[r['pos'] // W] = True
ph_div = np.zeros(nw, bool)
for p in P:
    se = os.path.join(P[p], 'phasing', f'{p}.switch_errors.bed')
    if os.path.exists(se):
        for line in open(se):
            f = line.split('\t')
            if len(f) > 1 and f[0] == a.chrom: ph_div[int(f[1]) // W] = True
bnd = {p: set((blocks_of(p)['from'].values // W).tolist()) for p in P}
for w in bnd['ont'] ^ bnd['hifi']: ph_div[w] = True
tab = [[int((sv_div & ph_div).sum()), int((sv_div & ~ph_div).sum())], [int((~sv_div & ph_div).sum()), int((~sv_div & ~ph_div).sum())]]
odds, pv = stats.fisher_exact(tab)
xa.append({'platform': 'both', 'test': 'overlap_of_SV_divergent_and_phasing_divergent_100kb_windows', 'n_gaps': int(sv_div.sum()), 'n_tp_sv': int(ph_div.sum()),
           'observed_fraction': tab[0][0] / max(sv_div.sum(), 1), 'null_mean_fraction': ph_div.mean(), 'null_sd': '', 'enrichment': odds, 'empirical_p': pv,
           'note': f'observed = fraction of SV-divergent windows that are also phasing-divergent; null = overall fraction of phasing-divergent windows; enrichment = Fisher odds ratio; table={tab}'})
pd.DataFrame(xa).to_csv(os.path.join(a.outdir, 'cross_arm.tsv'), sep='\t', index=False, float_format='%.4g')
print(pd.DataFrame(xa).to_string())
print(sz[sz.size_bin == 'ALL'][['svtype', 'ont_recall', 'hifi_recall', 'crossover']].to_string())
