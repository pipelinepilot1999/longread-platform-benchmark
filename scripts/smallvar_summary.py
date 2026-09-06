#!/usr/bin/env python3
"""Tidy hap.py extended.csv into SNP/INDEL x stratum rows (PASS-only calls)."""
import argparse, csv

ap = argparse.ArgumentParser()
ap.add_argument('--platform', required=True); ap.add_argument('--extended', required=True); ap.add_argument('--out', required=True)
a = ap.parse_args()
rows = []
for r in csv.DictReader(open(a.extended)):
    if r['Filter'] != 'PASS' or r['Subtype'] != '*' or r['Subset'] == '':
        continue
    rows.append({
        'platform': a.platform, 'type': r['Type'], 'stratum': r['Subset'],
        'truth_total': r['TRUTH.TOTAL'], 'truth_tp': r['TRUTH.TP'], 'truth_fn': r['TRUTH.FN'],
        'query_total': r['QUERY.TOTAL'], 'query_fp': r['QUERY.FP'], 'query_unk': r['QUERY.UNK'],
        'recall': r['METRIC.Recall'], 'precision': r['METRIC.Precision'], 'f1': r['METRIC.F1_Score'],
        'subset_size': r.get('Subset.Size', ''),
    })
with open(a.out, 'w') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter='\t'); w.writeheader(); w.writerows(rows)
for r in rows:
    if r['stratum'] == '*': print(r)
