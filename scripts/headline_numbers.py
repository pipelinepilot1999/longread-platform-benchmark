#!/usr/bin/env python3
"""Single source of truth for every number quoted in README.md and docs/.

Reads only files committed under results/ (standard library only, so it runs on any machine with
python3 and no network). Used by tests/test_docs_match_results.py and by `make numbers`.
"""
import csv, json, os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(ROOT, 'results')
BEST = {'ont': 'cutesv', 'hifi': 'sniffles'}
REGIMES = ['permissive', 'default', 'strict']
SIZE_BINS = ['50-99', '100-299', '300-999', '1k-9.9k', '10k-50k']
PRIMARY_DMRS = ['MCTS2P,HM13', 'BLCAP,NNAT', 'L3MBTL', 'GNAS', 'NESP-AS,GNAS-AS1', 'GNAS_XL', 'GNAS_Ex1A']   # the seven; the 5-CpG 'GNAS,GNAS-AS1' entry is reported separately


def tsv(*parts):
    with open(os.path.join(R, *parts)) as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def load():
    n = {}
    # ---- QC
    n['qc'] = {r['platform']: r for r in tsv('qc_summary.tsv')}
    # ---- SV benchmark
    sv = tsv('sv_benchmark.tsv')
    n['sv'] = {(r['platform'], r['caller'], r['regime'], r['svtype'], r['size_bin']): r for r in sv}
    n['sv_truth_in_bed'] = {t: int(n['sv'][('ont', 'cutesv', 'default', t, 'ALL')]['tp_base']) + int(n['sv'][('ont', 'cutesv', 'default', t, 'ALL')]['fn'])
                            for t in ('ALL', 'DEL', 'INS')}
    n['sv_concordance'] = {(r['platform'], r['category'], r['representative_caller'], r['svtype']): r for r in tsv('sv_concordance.tsv')}
    n['platform_concordance'] = {(r['category'], r['representative_platform'], r['svtype']): r for r in tsv('platform_concordance.tsv')}
    n['sv_by_size'] = {(r['svtype'], r['size_bin']): r for r in tsv('sv_platform_by_size.tsv')}
    with open(os.path.join(R, 'ont', 'truvari', 'platform_concordance', 'summary.json')) as fh:
        n['platform_truvari'] = json.load(fh)
    # ---- audit artefacts (2026-09-18): --pick ac sensitivity runs and recounts from the uncommitted VCFs
    n['sv_ac'] = {(r['platform'], r['caller'], r['regime']): r for r in tsv('sv_benchmark_pick_ac.tsv')}
    n['audit'] = {r['quantity']: r['value'] for r in tsv('audit_truth_counts.tsv')}
    n['cutesv_missed'] = tsv('audit_cutesv_hifi_missed_ins.tsv')
    # ---- small variants
    n['smallvar'] = {(r['platform'], r['type'], r['stratum']): r for r in tsv('smallvar_benchmark.tsv')}
    # ---- phasing
    n['phasing'] = {r['platform']: r for r in tsv('phasing_stats.tsv')}
    n['blocks'] = {p: tsv(p, 'phasing', f'{p}.blocks.tsv') for p in ('ont', 'hifi')}
    # ---- methylation (ONT)
    with open(os.path.join(R, 'ont', 'methylation_global.tsv')) as fh:
        head = fh.readline().rstrip('\n').split('\t'); vals = fh.readline().rstrip('\n').split('\t')
    n['meth_global'] = dict(zip(head, vals))
    n['meth_context'] = {(r['section'], r['category']): r for r in tsv('ont', 'methylation_by_context.tsv')}
    n['meth_sites'] = {int(r['threshold']): r for r in tsv('ont', 'methylation_site_counts.tsv')}
    ctrl = tsv('ont', 'imprinting_control.tsv')
    n['imprinting'] = {r['name']: r for r in ctrl if r['class'] != 'nanopore_detected'}
    n['imprinting_detected'] = {r['name']: r for r in ctrl if r['class'] == 'nanopore_detected'}
    n['imprinting_pass'] = sum(n['imprinting'][d]['control_pass'] == 'True' for d in PRIMARY_DMRS)
    def right_direction(r):
        # maternal DMRs of the GNAS cluster, NNAT and L3MBTL1 vs the paternal NESP55 DMR: direction is judged
        # relative to the locus' own neighbours, so here it just means a separation of at least 0.3.
        return float(r['haplotype_separation']) >= 0.3
    n['imprinting_directional'] = sum(right_direction(n['imprinting'][d]) for d in PRIMARY_DMRS)
    regions = tsv('ont', 'asm_regions.tsv')
    n['asm_regions'] = regions
    n['asm_sig_sites'] = sum(1 for _ in open(os.path.join(R, 'ont', 'asm_sites.tsv'))) - 1
    n['asm_tested'] = int(n['meth_sites'][5]['sites_both_haplotypes'])
    big = [r for r in regions if int(r['n_significant_cpg']) >= 20]
    n['asm_big_regions'] = len(big)
    n['asm_big_known'] = sum(r['classification'] == 'known_imprinted_control' for r in big)
    # ---- cross-arm
    n['cross'] = {(r['platform'], r['test']): r for r in tsv('cross_arm.tsv')}
    return n


def f1(n, platform, regime='default', caller=None, svtype='ALL', size='ALL', field='f1'):
    r = n['sv'][(platform, caller or BEST[platform], regime, svtype, size)]
    if field == 'recall':      # recompute from counts; the table stores 4 decimals and double rounding bites (0.6825 -> 0.683)
        return int(r['tp_base']) / (int(r['tp_base']) + int(r['fn']))
    if field == 'precision':
        return int(r['tp_comp']) / (int(r['tp_comp']) + int(r['fp']))
    return float(r[field])


def ac(n, platform, caller=None, regime='default', field='recall'):
    """--pick ac sensitivity run; rates recomputed from counts (the TSV's 4-decimal 0.7855 would double-round to 0.785)."""
    r = n['sv_ac'][(platform, caller or BEST[platform], regime)]
    rc = int(r['tp_base']) / (int(r['tp_base']) + int(r['fn'])); pr = int(r['tp_comp']) / (int(r['tp_comp']) + int(r['fp']))
    return {'recall': rc, 'precision': pr, 'f1': 2 * pr * rc / (pr + rc)}[field]


def headline(n=None):
    """The flat list of (label, formatted value) pairs that README, results.md and the guide must agree on."""
    n = n or load()
    q, p, s = n['qc'], n['phasing'], n['smallvar']
    out = [
        ('matched coverage target (x)', q['ont']['target_coverage']),
        ('native coverage ONT / HiFi (x)', f"{float(q['ont']['native_coverage']):.1f} / {float(q['hifi']['native_coverage']):.1f}"),
        ('effective chr20 depth ONT / HiFi (x)', f"{float(q['ont']['mean_depth_target_chrom']):.1f} / {float(q['hifi']['mean_depth_target_chrom']):.1f}"),
        ('read N50 ONT / HiFi (bp)', f"{int(q['ont']['read_n50']):,} / {int(q['hifi']['read_n50']):,}"),
        ('truth SVs in benchmark BED (all / DEL / INS)', ' / '.join(str(n['sv_truth_in_bed'][t]) for t in ('ALL', 'DEL', 'INS'))),
        ('SV F1 default, ONT cuteSV / HiFi Sniffles2', f"{f1(n,'ont'):.3f} / {f1(n,'hifi'):.3f}"),
        ('SV recall default, ONT / HiFi', f"{f1(n,'ont',field='recall'):.3f} / {f1(n,'hifi',field='recall'):.3f}"),
        ('SV recall default with --pick ac, ONT / HiFi', f"{ac(n,'ont'):.3f} / {ac(n,'hifi'):.3f}"),
        ('SV F1 default with --pick ac, ONT / HiFi', f"{ac(n,'ont',field='f1'):.3f} / {ac(n,'hifi',field='f1'):.3f}"),
        ('SV precision default, ONT / HiFi', f"{f1(n,'ont',field='precision'):.3f} / {f1(n,'hifi',field='precision'):.3f}"),
        ('DEL 50-99 bp recall ONT / HiFi', f"{f1(n,'ont',svtype='DEL',size='50-99',field='recall'):.3f} / {f1(n,'hifi',svtype='DEL',size='50-99',field='recall'):.3f}"),
        ('INS recall default -> strict, ONT', f"{f1(n,'ont',svtype='INS',field='recall'):.3f} -> {f1(n,'ont','strict',svtype='INS',field='recall'):.3f}"),
        ('INS recall default -> strict, HiFi', f"{f1(n,'hifi',svtype='INS',field='recall'):.3f} -> {f1(n,'hifi','strict',svtype='INS',field='recall'):.3f}"),
        ('chr20 truth SVs 50 bp-50 kb (all / DEL / INS)', f"{int(n['audit']['truth_sv_chr20_50bp_50kb_all']):,} / {n['audit']['truth_sv_chr20_50bp_50kb_DEL']} / {int(n['audit']['truth_sv_chr20_50bp_50kb_INS']):,}"),
        ('truth hets within 45 kb of L3MBTL1 DMR / nearest', f"{n['audit']['truth_het_within_45kb_of_L3MBTL1_DMR']} / {int(n['audit']['nearest_truth_het_left_of_L3MBTL1_DMR_distance_bp'])/1000:.1f} kb"),
        ('SNV F1 ONT / HiFi', f"{float(s[('ont','SNP','*')]['f1']):.4f} / {float(s[('hifi','SNP','*')]['f1']):.4f}"),
        ('indel F1 ONT / HiFi', f"{float(s[('ont','INDEL','*')]['f1']):.3f} / {float(s[('hifi','INDEL','*')]['f1']):.3f}"),
        ('indel F1 in homopolymers >= 12 bp ONT / HiFi', f"{float(s[('ont','INDEL','homopolymer_ge12')]['f1']):.3f} / {float(s[('hifi','INDEL','homopolymer_ge12')]['f1']):.3f}"),
        ('phase block N50 ONT / HiFi (Mb)', f"{int(p['ont']['block_n50_bp'])/1e6:.2f} / {int(p['hifi']['block_n50_bp'])/1e6:.2f}"),
        ('switch error ONT / HiFi (%)', f"{100*float(p['ont']['switch_error_rate']):.3f} / {100*float(p['hifi']['switch_error_rate']):.3f}"),
        ('reads unassigned ONT / HiFi (%)', f"{100*float(p['ont']['fraction_reads_unassigned']):.1f} / {100*float(p['hifi']['fraction_reads_unassigned']):.1f}"),
        ('CpGs >= 10x unphased', f"{int(n['meth_global']['n_sites']):,}"),
        ('promoter CpG island mean methylation', f"{float(n['meth_context'][('promoter_islands','promoter_cpg_island')]['mean_methylation']):.3f}"),
        ('imprinting control: primary DMRs passing control_pass', f"{n['imprinting_pass']} of {len(PRIMARY_DMRS)}"),
        ('ASM: significant / tested CpGs', f"{n['asm_sig_sites']:,} / {n['asm_tested']:,}"),
        ('ASM regions (>= 5 CpGs) / with >= 20 CpGs / of those known imprinted', f"{len(n['asm_regions']):,} / {n['asm_big_regions']} / {n['asm_big_known']}"),
    ]
    return out


if __name__ == '__main__':
    for k, v in headline():
        print(f'{k:68s} {v}')
