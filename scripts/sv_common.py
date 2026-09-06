"""Shared helpers for SV tables: size bins, SV type/length extraction from VCF records."""
import pysam

SIZE_BINS = [(50, 100, '50-99'), (100, 300, '100-299'), (300, 1000, '300-999'), (1000, 10000, '1k-9.9k'), (10000, 50001, '10k-50k')]

def size_bin(l):
    for lo, hi, name in SIZE_BINS:
        if lo <= l < hi:
            return name
    return 'other'

def sv_len(rec):
    info = rec.info
    if 'SVLEN' in info:
        v = info['SVLEN']
        v = v[0] if isinstance(v, tuple) else v
        return abs(int(v))
    alt = rec.alts[0] if rec.alts else ''
    if alt.startswith('<'):
        return abs(rec.stop - rec.start)
    return abs(len(alt) - len(rec.ref))

def sv_type(rec):
    t = rec.info.get('SVTYPE')
    if t:
        return t if isinstance(t, str) else t[0]
    alt = rec.alts[0] if rec.alts else ''
    return 'INS' if len(alt) > len(rec.ref) else 'DEL'

def read_vcf(path, extra=()):
    rows = []
    with pysam.VariantFile(path) as vf:
        for rec in vf:
            d = {'chrom': rec.chrom, 'pos': rec.pos, 'id': rec.id, 'svtype': sv_type(rec), 'svlen': sv_len(rec)}
            d['size_bin'] = size_bin(d['svlen'])
            for k in extra:
                v = rec.info.get(k)
                d[k] = v[0] if isinstance(v, tuple) else v
            rows.append(d)
    return rows
