#!/usr/bin/env python3
"""Build a hap.py stratification TSV for one chromosome from the GIAB v3.6 stratifications.

Homopolymer strata are the point of this: ONT error concentrates in homopolymers, so indel
performance is reported inside and outside them. Every BED is subset to the chromosome so
hap.py does not spend time on the rest of the genome.
"""
import argparse, gzip, os

STRATA = {
    # name: relative path under the GIAB GRCh38@all directory
    'homopolymer_4to6':        'LowComplexity/GRCh38_SimpleRepeat_homopolymer_4to6_slop5.bed.gz',
    'homopolymer_7to11':       'LowComplexity/GRCh38_SimpleRepeat_homopolymer_7to11_slop5.bed.gz',
    'homopolymer_ge12':        'LowComplexity/GRCh38_SimpleRepeat_homopolymer_ge12_slop5.bed.gz',
    'imperfect_homopol_ge11':  'LowComplexity/GRCh38_SimpleRepeat_imperfecthomopolge11_slop5.bed.gz',
    'all_homopolymers_ge7':    'LowComplexity/GRCh38_AllHomopolymers_ge7bp_imperfectge11bp_slop5.bed.gz',
    'not_in_homopolymers_ge7': 'LowComplexity/GRCh38_notinAllHomopolymers_ge7bp_imperfectge11bp_slop5.bed.gz',
    'all_TR_and_homopolymers': 'LowComplexity/GRCh38_AllTandemRepeatsandHomopolymers_slop5.bed.gz',
    'not_in_TR_or_homopolymers': 'LowComplexity/GRCh38_notinAllTandemRepeatsandHomopolymers_slop5.bed.gz',
    'segdups':                 'SegmentalDuplications/GRCh38_segdups.bed.gz',
    'not_in_segdups':          'SegmentalDuplications/GRCh38_notinsegdups.bed.gz',
}

ap = argparse.ArgumentParser()
ap.add_argument('--giab', required=True); ap.add_argument('--chrom', required=True); ap.add_argument('--outdir', default='.')
a = ap.parse_args()
rows = []
for name, rel in STRATA.items():
    src = os.path.join(a.giab, rel)
    if not os.path.exists(src):
        print(f'WARNING missing {src}; stratum {name} skipped'); continue
    out = os.path.join(a.outdir, f'{name}.{a.chrom}.bed.gz')
    n = 0
    with gzip.open(src, 'rt') as fh, gzip.open(out, 'wt') as oh:
        for line in fh:
            if line.startswith(a.chrom + '\t'):
                oh.write(line); n += 1
    rows.append((name, os.path.basename(out), n))
with open(os.path.join(a.outdir, 'strata.tsv'), 'w') as fh:
    for name, f, n in rows:
        fh.write(f'{name}\t{f}\n')
for r in rows: print(*r)
