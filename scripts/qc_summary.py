#!/usr/bin/env python3
"""Collate per-platform QC into one row: reads, bases, N50, native/matched coverage, mapping and
supplementary rates, fraction of pre-selected reads that realigned to the target chromosome, mean
depth, and MM/ML tag survival."""
import argparse, csv, re

ap = argparse.ArgumentParser()
for k in ['platform', 'native', 'downsample', 'seqkit', 'flagstat', 'idxstats', 'mosdepth', 'tags', 'chrom', 'out']:
    ap.add_argument('--' + k.replace('_', '-'), required=True)
a = ap.parse_args()

nat = next(csv.DictReader(open(a.native), delimiter='\t'))
ds = next(csv.DictReader(open(a.downsample), delimiter='\t'))
sk = next(csv.DictReader(open(a.seqkit), delimiter='\t'))
tags = next(csv.DictReader(open(a.tags), delimiter='\t'))

fs = {}
for line in open(a.flagstat):
    m = re.match(r'(\d+) \+ \d+ (.+)', line.strip())
    if m:
        key = m.group(2).split(' (')[0]
        fs[key] = int(m.group(1))
total = fs.get('in total', 0); secondary = fs.get('secondary', 0); supp = fs.get('supplementary', 0)
mapped = fs.get('mapped', 0)
primary_total = total - secondary - supp
primary_mapped = mapped - secondary - supp

idx = {r[0]: int(r[2]) for r in csv.reader(open(a.idxstats), delimiter='\t')}
on_chrom = idx.get(a.chrom, 0); all_mapped = sum(idx.values())

md = {}
for r in csv.DictReader(open(a.mosdepth), delimiter='\t'):
    md[r['chrom']] = r
mean_depth = float(md[a.chrom]['mean']) if a.chrom in md else float('nan')

row = {
    'platform': a.platform,
    'raw_reads': nat['reads'], 'raw_bases': nat['bases'], 'native_coverage': nat['native_coverage'],
    'downsample_mode': ds['mode'], 'target_coverage': ds['target_coverage'], 'downsample_fraction': ds['fraction'],
    'reads_analysed': sk['num_seqs'], 'bases_analysed': sk['sum_len'], 'read_n50': sk['N50'],
    'read_mean': sk['avg_len'], 'read_max': sk['max_len'], 'read_min': sk['min_len'],
    'primary_records': primary_total, 'primary_mapped': primary_mapped,
    'mapping_rate': f'{primary_mapped/primary_total:.4f}' if primary_total else 'NA',
    'supplementary_records': supp, 'supplementary_per_primary': f'{supp/primary_total:.4f}' if primary_total else 'NA',
    'records_on_target_chrom': on_chrom, 'fraction_records_on_target_chrom': f'{on_chrom/all_mapped:.4f}' if all_mapped else 'NA',
    'mean_depth_target_chrom': f'{mean_depth:.2f}',
    'tagged_MM_fraction': tags['fraction_MM'], 'tagged_MM_reads': tags['with_MM'],
}
with open(a.out, 'w') as fh:
    fh.write('\t'.join(row.keys()) + '\n'); fh.write('\t'.join(str(v) for v in row.values()) + '\n')
print(row)
