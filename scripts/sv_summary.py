#!/usr/bin/env python3
"""Per-platform SV summary from Truvari outputs.

sv_benchmark.<platform>.tsv   caller x regime x svtype x size bin: TP/FP/FN, precision, recall, F1
sv_concordance.<platform>.tsv shared / caller-unique calls (inside the benchmark BED) cross-referenced
                              against their default-regime truth status; answers whether caller-unique
                              calls are enriched for false positives.
sv_calls_annotated.<platform>.tsv.gz  every benchmarked call with its status under each regime.
"""
import argparse, csv, gzip, json, os, sys
from collections import Counter, defaultdict
sys.path.insert(0, os.path.dirname(__file__))
from sv_common import read_vcf, SIZE_BINS

ap = argparse.ArgumentParser()
ap.add_argument('--platform', required=True); ap.add_argument('--bench', nargs='+', required=True)
ap.add_argument('--concordance', required=True); ap.add_argument('--outdir', default='.')
a = ap.parse_args()

def prf(tp, fp, fn):
    p = tp / (tp + fp) if tp + fp else 0.0
    r = tp / (tp + fn) if tp + fn else 0.0
    f = 2 * p * r / (p + r) if p + r else 0.0
    return p, r, f

bench_rows = []
status = defaultdict(dict)          # (caller, id) -> {regime: TP/FP}
call_info = {}
for d in a.bench:
    caller, regime = os.path.basename(d).split('_', 1)
    summ = json.load(open(os.path.join(d, 'summary.json')))
    tpb = read_vcf(os.path.join(d, 'tp-base.vcf.gz')); fn = read_vcf(os.path.join(d, 'fn.vcf.gz'))
    tpc = read_vcf(os.path.join(d, 'tp-comp.vcf.gz')); fp = read_vcf(os.path.join(d, 'fp.vcf.gz'))
    for r in tpc: status[(caller, r['id'])][regime] = 'TP'; call_info[(caller, r['id'])] = r
    for r in fp:  status[(caller, r['id'])][regime] = 'FP'; call_info[(caller, r['id'])] = r
    # overall + stratified
    strata = [('ALL', 'ALL')] + [(t, 'ALL') for t in ('DEL', 'INS')] + [(t, b[2]) for t in ('DEL', 'INS') for b in SIZE_BINS]
    for t, b in strata:
        sel = lambda rows: [r for r in rows if (t == 'ALL' or r['svtype'] == t) and (b == 'ALL' or r['size_bin'] == b)]
        ntp_b, nfn, ntp_c, nfp = len(sel(tpb)), len(sel(fn)), len(sel(tpc)), len(sel(fp))
        p, r, f = prf(ntp_c, nfp, nfn)
        recall = ntp_b / (ntp_b + nfn) if ntp_b + nfn else 0.0
        f = 2 * p * recall / (p + recall) if p + recall else 0.0
        bench_rows.append({'platform': a.platform, 'caller': caller, 'regime': regime, 'svtype': t, 'size_bin': b,
                           'tp_base': ntp_b, 'fn': nfn, 'tp_comp': ntp_c, 'fp': nfp,
                           'precision': f'{p:.4f}', 'recall': f'{recall:.4f}', 'f1': f'{f:.4f}',
                           'truvari_precision': f"{summ['precision']:.4f}" if (t, b) == ('ALL', 'ALL') else '',
                           'truvari_recall': f"{summ['recall']:.4f}" if (t, b) == ('ALL', 'ALL') else '',
                           'truvari_f1': f"{summ['f1']:.4f}" if (t, b) == ('ALL', 'ALL') else ''})
bench_rows.sort(key=lambda r: (r['caller'], ['permissive', 'default', 'strict'].index(r['regime']) if r['regime'] in ('permissive', 'default', 'strict') else 9, r['svtype'], r['size_bin']))
with open(os.path.join(a.outdir, f'sv_benchmark.{a.platform}.tsv'), 'w') as fh:
    w = csv.DictWriter(fh, fieldnames=list(bench_rows[0].keys()), delimiter='\t'); w.writeheader(); w.writerows(bench_rows)

# ---- caller concordance (base = sniffles, comp = cutesv, default regime, inside the benchmark BED)
cd = a.concordance
shared_s = read_vcf(os.path.join(cd, 'tp-base.vcf.gz'), extra=('MatchId',))
only_s   = read_vcf(os.path.join(cd, 'fn.vcf.gz'))
shared_c = read_vcf(os.path.join(cd, 'tp-comp.vcf.gz'), extra=('MatchId',))
only_c   = read_vcf(os.path.join(cd, 'fp.vcf.gz'))
def truth(caller, rid):
    return status.get((caller, rid), {}).get('default', 'not_benchmarked')
conc_rows = []
groups = [('shared', 'sniffles', shared_s), ('shared', 'cutesv', shared_c), ('sniffles_only', 'sniffles', only_s), ('cutesv_only', 'cutesv', only_c)]
for cat, caller, rows in groups:
    for t in ('ALL', 'DEL', 'INS'):
        sel = [r for r in rows if t == 'ALL' or r['svtype'] == t]
        c = Counter(truth(caller, r['id']) for r in sel)
        n = len(sel); tp = c.get('TP', 0); fp = c.get('FP', 0)
        conc_rows.append({'platform': a.platform, 'category': cat, 'representative_caller': caller, 'svtype': t, 'n_calls': n,
                          'truth_TP': tp, 'truth_FP': fp, 'not_benchmarked': c.get('not_benchmarked', 0),
                          'fp_fraction_of_benchmarked': f'{fp/(tp+fp):.4f}' if tp + fp else 'NA'})
with open(os.path.join(a.outdir, f'sv_concordance.{a.platform}.tsv'), 'w') as fh:
    w = csv.DictWriter(fh, fieldnames=list(conc_rows[0].keys()), delimiter='\t'); w.writeheader(); w.writerows(conc_rows)

shared_ids = {('sniffles', r['id']) for r in shared_s} | {('cutesv', r['id']) for r in shared_c}
with gzip.open(os.path.join(a.outdir, f'sv_calls_annotated.{a.platform}.tsv.gz'), 'wt') as fh:
    fh.write('platform\tcaller\tid\tchrom\tpos\tsvtype\tsvlen\tsize_bin\tshared_between_callers\tstatus_permissive\tstatus_default\tstatus_strict\n')
    for (caller, rid), r in sorted(call_info.items(), key=lambda kv: (kv[0][0], kv[1]['pos'])):
        s = status[(caller, rid)]
        fh.write(f"{a.platform}\t{caller}\t{rid}\t{r['chrom']}\t{r['pos']}\t{r['svtype']}\t{r['svlen']}\t{r['size_bin']}\t{(caller, rid) in shared_ids}\t{s.get('permissive','')}\t{s.get('default','')}\t{s.get('strict','')}\n")
for r in bench_rows:
    if r['svtype'] == 'ALL': print(r['caller'], r['regime'], 'P', r['precision'], 'R', r['recall'], 'F1', r['f1'])
for r in conc_rows:
    if r['svtype'] == 'ALL': print(r)
