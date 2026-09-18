#!/usr/bin/env python3
"""Audit recounts that need the large, uncommitted files (truth VCFs, call VCFs, Truvari outputs).

Run on the host inside the analysis image (pysam):
  docker run --rm -u $(id -u):$(id -g) -v /home/ec2-user:/home/ec2-user longread-platform-benchmark/analysis:1.0 \
      python3 /home/ec2-user/longread/repo/scripts/audit_recounts.py
Writes results/audit_truth_counts.tsv (truth SV totals, heterozygous truth records, heterozygosity around the
L3MBTL1 DMR) and results/audit_cutesv_hifi_missed_ins.tsv (why each truth insertion >= 1 kb that cuteSV 'missed'
on HiFi was scored as a miss)."""
import pysam, sys, os, bisect, csv
R = '/home/ec2-user/longread/resources'; REPO = '/home/ec2-user/longread/repo'
CH = 'chr20'; rows = []
def add(k, v, note=''): rows.append((k, v, note)); print(k, v, note, sep='\t')

def svinfo(rec):
    alt = rec.alts[0] if rec.alts else ''
    L = rec.info.get('SVLEN'); L = L[0] if isinstance(L, tuple) else L
    if L is None: L = len(alt) - len(rec.ref)
    t = rec.info.get('SVTYPE'); t = t if isinstance(t, str) or t is None else t[0]
    if not t: t = 'INS' if len(alt) > len(rec.ref) else 'DEL'
    return t, abs(int(L))

bed = [(int(l.split()[1]), int(l.split()[2])) for l in open(f'{R}/truth_chr20/truth_sv.chr20.bed')]
bed.sort(); starts = [b[0] for b in bed]
def within_bed(s, e):          # truvari --includebed: variant fully contained in one region
    i = bisect.bisect_right(starts, s) - 1
    return i >= 0 and bed[i][0] <= s and e <= bed[i][1]
add('sv_bed_regions_chr20', len(bed)); add('sv_bed_bp_chr20', sum(e - s for s, e in bed))

tot = {'all': 0, 'DEL': 0, 'INS': 0}; inb = {'all': 0, 'DEL': 0, 'INS': 0}; ge50 = 0; gt50k = 0
for rec in pysam.VariantFile(f'{R}/truth_chr20/truth_sv.chr20.vcf.gz'):
    t, L = svinfo(rec)
    if L < 50: continue
    ge50 += 1
    if L > 50000: gt50k += 1; continue
    tot['all'] += 1; tot[t] = tot.get(t, 0) + 1
    if within_bed(rec.start, rec.stop): inb['all'] += 1; inb[t] = inb.get(t, 0) + 1
add('truth_sv_chr20_ge50bp_after_ALTstar_removal', ge50); add('truth_sv_chr20_gt50kb', gt50k)
for k in ('all', 'DEL', 'INS'): add(f'truth_sv_chr20_50bp_50kb_{k}', tot[k])
for k in ('all', 'DEL', 'INS'): add(f'truth_sv_chr20_50bp_50kb_inside_bed_{k}', inb[k], 'python containment check; Truvari reports 802/254/548')

# ---- small-variant truth (FILTER '.' only, as used by the pipeline)
n = het = 0; hetpos = []
for rec in pysam.VariantFile(f'{R}/truth_chr20/truth_sm.chr20.vcf.gz'):
    n += 1; gt = rec.samples[0]['GT']
    if gt and None not in gt and len(set(gt)) > 1: het += 1; hetpos.append(rec.pos)
add('truth_smvar_chr20_records_filter_dot', n); add('truth_smvar_chr20_heterozygous', het)
add('truth_het_per_Mb_chr20', round(het / 64.444167, 1))
hetpos.sort()
def count(a, lo, hi): return bisect.bisect_right(a, hi) - bisect.bisect_left(a, lo)
D0, D1 = 43513725, 43515400      # L3MBTL1 DMR
for pad in (10000, 20000, 45000, 100000, 500000):
    add(f'truth_het_within_{pad//1000}kb_of_L3MBTL1_DMR', count(hetpos, D0 - pad, D1 + pad))
left = [p for p in hetpos if p < D0]; right = [p for p in hetpos if p > D1]
add('nearest_truth_het_left_of_L3MBTL1_DMR_distance_bp', D0 - left[-1] if left else 'NA'); add('nearest_truth_het_right_distance_bp', right[0] - D1 if right else 'NA')
# unfiltered v5.0q (incl. HET/GAP flagged records) as a cross-check
allhet = []
for rec in pysam.VariantFile(f'{R}/v5.0q/HG002_GRCh38_v5.0q_smvar.vcf.gz').fetch(CH, D0 - 500000, D1 + 500000):
    gt = rec.samples[0]['GT']
    if gt and None not in gt and len(set(gt)) > 1: allhet.append(rec.pos)
add('unfiltered_v5.0q_het_within_45kb_of_L3MBTL1_DMR', count(sorted(allhet), D0 - 45000, D1 + 45000), 'includes HET/GAP-flagged records')
for p in ('ont', 'hifi'):
    pos = []
    for rec in pysam.VariantFile(f'{REPO}/results/{p}/phasing/{p}.phased.vcf.gz').fetch(CH, D0 - 45000, D1 + 45000):
        s = rec.samples[0]; gt = s['GT']
        if gt and None not in gt and len(set(gt)) > 1: pos.append((rec.pos, s.phased, rec.qual, rec.ref, rec.alts[0][:12]))
    add(f'{p}_clair3_pass_het_within_45kb_of_L3MBTL1_DMR', len(pos), 'phased=%d; ' % sum(1 for x in pos if x[1]) + ';'.join(f'{x[0]}:{x[3][:6]}>{x[4]}:q{x[2]:.0f}' for x in pos[:12]))

with open(f'{REPO}/results/audit_truth_counts.tsv', 'w') as fh:
    fh.write('quantity\tvalue\tnote\n'); [fh.write(f'{a}\t{b}\t{c}\n') for a, b, c in rows]

# ================= cuteSV on HiFi: why is each truth insertion >= 1 kb scored as missed?
T = f'{REPO}/results/hifi/truvari'
def L(rec):
    return svinfo(rec)[1]
def load(path): return list(pysam.VariantFile(path))
def is_ins(r): return bool(r.alts) and (len(r.alts[0]) > len(r.ref) or r.alts[0].startswith('<INS'))
trs_all = sorted((int(l.split()[1]), int(l.split()[2])) for l in open(f'{R}/human_GRCh38_no_alt_analysis_set.trf.bed') if l.startswith(CH + '\t'))
tr_starts = [x[0] for x in trs_all]
def in_tr(pos, pad=50):
    i = bisect.bisect_right(tr_starts, pos + pad) - 1
    while i >= 0 and trs_all[i][0] > pos - 100000:
        if trs_all[i][0] - pad <= pos <= trs_all[i][1] + pad: return True
        i -= 1
    return False
fn = [r for r in load(f'{T}/cutesv_default/fn.vcf.gz') if L(r) >= 1000 and is_ins(r)]
tpb = load(f'{T}/cutesv_default/tp-base.vcf.gz'); tpc = load(f'{T}/cutesv_default/tp-comp.vcf.gz'); fpc = load(f'{T}/cutesv_default/fp.vcf.gz')
perm = {(r.pos, L(r)) for r in load(f'{T}/cutesv_permissive/tp-base.vcf.gz')}
ac = {(r.pos, L(r)) for r in load(f'{REPO}/results/hifi/truvari_pick_ac/cutesv_default/tp-base.vcf.gz')} if os.path.exists(f'{REPO}/results/hifi/truvari_pick_ac/cutesv_default/tp-base.vcf.gz') else set()
snif_fn = {(r.pos, L(r)) for r in load(f'{T}/sniffles_default/fn.vcf.gz')}
out = []
for r in fn:
    tl = L(r)
    near_tp = [c for c in tpc if abs(c.pos - r.pos) <= 1000 and is_ins(c)]
    sister = [b for b in tpb if abs(b.pos - r.pos) <= 1000 and L(b) >= 50]
    near_fp = [c for c in fpc if abs(c.pos - r.pos) <= 1000]
    if near_tp and sister: why = 'two_truth_records_one_call'      # the single cuteSV call was consumed by the other truth record (--pick single)
    elif near_tp: why = 'call_matched_other_truth_record'
    elif near_fp: why = 'call_present_but_scored_FP'
    else: why = 'no_cutesv_call_within_1kb'
    g = lambda k: (round(float(r.info[k]), 4) if k in r.info else '')
    out.append(dict(pos=r.pos, truth_len=tl, truth_gt='/'.join(map(str, r.samples[0]['GT'])), in_tandem_repeat=in_tr(r.pos), why_missed=why,
                    best_candidate_seq_similarity=g('PctSeqSimilarity'), best_candidate_size_similarity=g('PctSizeSimilarity'),
                    recovered_under_permissive=(r.pos, tl) in perm, recovered_with_pick_ac=(r.pos, tl) in ac, also_missed_by_sniffles=(r.pos, tl) in snif_fn))
with open(f'{REPO}/results/audit_cutesv_hifi_missed_ins.tsv', 'w') as fh:
    w = csv.DictWriter(fh, fieldnames=list(out[0].keys()), delimiter='\t'); w.writeheader(); w.writerows(out)
from collections import Counter
print('\ncuteSV HiFi, truth INS >= 1 kb scored as missed:', len(out), Counter(o['why_missed'] for o in out))
print('in tandem repeat:', sum(o['in_tandem_repeat'] for o in out), '| recovered with --pick ac:', sum(o['recovered_with_pick_ac'] for o in out),
      '| recovered under permissive:', sum(o['recovered_under_permissive'] for o in out), '| also missed by Sniffles2:', sum(o['also_missed_by_sniffles'] for o in out))
