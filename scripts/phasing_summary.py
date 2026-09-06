#!/usr/bin/env python3
"""Phase-block statistics, fraction of heterozygous variants phased, switch/Hamming error versus the
assembly-phased truth, haplotag rates, and the phase-block-N50 vs read-N50 ratio (the section 12 assertion)."""
import argparse, csv, gzip, numpy as np

ap = argparse.ArgumentParser()
for k in ['platform', 'stats', 'blocks', 'compare', 'haplotag', 'read_lengths', 'out']:
    ap.add_argument('--' + k.replace('_', '-'), required=True)
ap.add_argument('--chrom-len', type=int, default=64444167, help='for the chromosome-relative N50 (whatshap convention)')
a = ap.parse_args()

st = [r for r in csv.DictReader(open(a.stats), delimiter='\t')]
st = [r for r in st if r['chromosome'] == 'ALL'][0] if any(r['chromosome'] == 'ALL' for r in st) else st[0]
blocks = [r for r in csv.DictReader(open(a.blocks), delimiter='\t')]
lengths = np.array([int(r['to']) - int(r['from']) + 1 for r in blocks]) if blocks else np.array([0])
def n50(x, total=None):
    x = np.sort(x)[::-1]; c = np.cumsum(x); total = total or c[-1]
    i = np.searchsorted(c, total / 2)
    return int(x[i]) if len(x) and i < len(x) else 0
rl = np.array([int(l) for l in gzip.open(a.read_lengths, 'rt') if l.strip()])
cmp = next(csv.DictReader(open(a.compare), delimiter='\t'))
ht = next(csv.DictReader(open(a.haplotag), delimiter='\t'))

het = int(st['heterozygous_variants']); phased = int(st['phased'])
row = {
    'platform': a.platform,
    'heterozygous_variants': het, 'phased_variants': phased, 'fraction_het_phased': f'{phased/het:.4f}' if het else 'NA',
    'phase_blocks': len(blocks), 'block_n50_bp': n50(lengths, a.chrom_len), 'block_n50_phased_span_bp': n50(lengths), 'phased_span_bp': int(lengths.sum()), 'fraction_chrom_in_blocks': f'{lengths.sum()/a.chrom_len:.4f}', 'block_max_bp': int(lengths.max()),
    'block_median_bp': int(np.median(lengths)), 'whatshap_block_n50_bp': st.get('block_n50', ''),
    'read_n50_bp': n50(rl), 'block_n50_over_read_n50': f'{n50(lengths, a.chrom_len)/n50(rl):.2f}' if n50(rl) else 'NA',
    'compare_common_hets': cmp.get('all_assessed_pairs', ''), 'switch_errors': cmp.get('all_switches', ''),
    'switch_error_rate': cmp.get('all_switch_rate', ''), 'switch_flips': cmp.get('all_switchflips', ''),
    'hamming_rate': cmp.get('blockwise_hamming_rate', ''), 'covered_variants_truth_and_query': cmp.get('covered_variants', ''),
    'primary_reads': ht['primary_reads'], 'reads_hp1': ht['hp1'], 'reads_hp2': ht['hp2'], 'reads_unassigned': ht['unassigned'],
    'fraction_reads_tagged': ht['fraction_tagged'], 'fraction_reads_unassigned': f"{int(ht['unassigned'])/int(ht['primary_reads']):.4f}",
}
with open(a.out, 'w') as fh:
    fh.write('\t'.join(row.keys()) + '\n'); fh.write('\t'.join(str(v) for v in row.values()) + '\n')
print(row)
