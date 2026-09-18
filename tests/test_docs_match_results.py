"""Docs-vs-results consistency tests (standard library + pytest/unittest; no network, no data files).

Every test here exists because the 2026-09-18 audit found the corresponding defect:
  - README shipped with unfilled {{PLACEHOLDER}} tokens
  - the imprinting control was reported as 5/7 while results/ont/imprinting_control.tsv says 4/7
  - prose claimed strict matching costs both platforms the same insertion recall (it does not)
  - decision ids (D-0xx) cited in code and docs pointed at the wrong or a non-existent entry
  - docs cited a modkit flag the code does not use, and file paths that do not exist
  - README commands used a --target_cov that differs from the one recorded in results/
"""
import glob, os, re, sys, unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'scripts'))
import headline_numbers as hn   # noqa: E402

N = hn.load()
DOCS = ['README.md', 'docs/results.md', 'docs/decisions.md', 'docs/limitations.md', 'docs/validation.md']


def read(p):
    with open(os.path.join(ROOT, p)) as fh:
        return fh.read()


def squash(s):
    return re.sub(r'\s+', ' ', s)


class Placeholders(unittest.TestCase):
    def test_no_unfilled_template_tokens(self):
        for d in DOCS:
            self.assertEqual(re.findall(r'\{\{[A-Z_]+\}\}', read(d)), [], f'unfilled placeholder in {d}')

    def test_no_pending_markers(self):
        for d in DOCS:
            self.assertNotIn('*pending*', read(d), d)


class Paths(unittest.TestCase):
    def test_backticked_repo_paths_exist(self):
        pat = re.compile(r'`((?:results|scripts|modules|docs|containers|tests)/[^`\s]+)`')
        missing = []
        for d in DOCS:
            for m in pat.findall(read(d)):
                path = m.rstrip('.,;:)')
                if any(tok in path for tok in ('<', '{', '*.')) and not glob.glob(os.path.join(ROOT, path.replace('{ont,hifi}', '*'))):
                    if '<' in path or '{' in path:
                        continue
                if path.startswith('results/') and path.endswith(('.bam', '.vcf.gz', '.bed.gz', '.tsv.gz')):
                    continue        # large outputs are deliberately not committed
                if not glob.glob(os.path.join(ROOT, path)) and not glob.glob(os.path.join(ROOT, path + '*')):
                    missing.append((d, path))
        self.assertEqual(missing, [])

    def test_readme_listed_scripts_exist(self):
        block = read('README.md').split('scripts/   ')[1].split('results/')[0]
        for name in re.findall(r'[\w]+\.(?:py|sh)', block):
            self.assertTrue(os.path.exists(os.path.join(ROOT, 'scripts', name)), name)


class DecisionIds(unittest.TestCase):
    def setUp(self):
        self.entries = dict(re.findall(r'^### (D-\d{3}): (.+)$', read('docs/decisions.md'), re.M))

    def test_every_cited_decision_exists(self):
        for f in DOCS + ['main.nf', 'nextflow.config'] + [os.path.relpath(p, ROOT) for p in glob.glob(os.path.join(ROOT, 'modules', '*.nf'))]:
            for d in set(re.findall(r'D-\d{3}', read(f))):
                self.assertIn(d, self.entries, f'{f} cites {d}, which is not in decisions.md')

    def test_topical_citations_point_at_the_right_entry(self):
        self.assertIn('One platform per run', self.entries['D-007'])
        self.assertIn('D-007', read('main.nf').split('nextflow.enable.dsl')[0])
        self.assertIn('Methylation on ONT only', self.entries['D-023'])
        self.assertRegex(read('README.md'), r'Methylation was run on ONT only[^\n]*\(D-023\)')


class CodeVsDocs(unittest.TestCase):
    def test_modkit_flags_in_docs_are_the_flags_in_code(self):
        code = read('modules/modkit_pileup.nf')
        self.assertIn('--modified-bases 5mC', code)
        self.assertNotIn('--ignore h', code)
        self.assertNotIn('--ignore h', read('docs/limitations.md'))

    def test_truvari_regimes_in_docs_match_config(self):
        cfg = read('nextflow.config')
        for regime in ("--refdist 1000 --pctseq 0 --pctsize 0.5 --typeignore", "--refdist 500 --pctseq 0.7 --pctsize 0.7", "--refdist 100 --pctseq 0.9 --pctsize 0.9"):
            self.assertIn(regime, cfg); self.assertIn(regime, read('docs/results.md')); self.assertIn(regime, read('docs/decisions.md'))

    def test_recorded_commands_match_config(self):
        self.assertIn('-ax map-ont -y', read('results/ont/qc/ont.minimap2.cmd'))
        self.assertIn('-ax map-hifi -y', read('results/hifi/qc/hifi.minimap2.cmd'))
        cfg = read('nextflow.config')
        for p in ('ont', 'hifi'):
            args = re.search(p + r"\s*:\s*'(--max_cluster_bias_INS[^']+)'", cfg).group(1)
            self.assertIn(args, read(f'results/{p}/sv/{p}.cutesv.cmd'))

    def test_readme_cli_params_exist_in_config(self):
        cfg = read('nextflow.config')
        for flag in set(re.findall(r'--([a-z_]+)', read('README.md').split('## Pipeline')[1].split('```')[1])):
            self.assertRegex(cfg, r'\b' + flag + r'\s*=', f'README uses --{flag}, not a param in nextflow.config')

    def test_readme_target_cov_is_the_one_that_was_run(self):
        self.assertIn(f"--target_cov {N['qc']['ont']['target_coverage']}", read('README.md'))

    def test_no_dead_params(self):
        cfg = read('nextflow.config')
        body = ''.join(read(os.path.relpath(p, ROOT)) for p in glob.glob(os.path.join(ROOT, 'modules', '*.nf'))) + read('main.nf')
        block = cfg.split('params {')[1].split('\nprocess {')[0]
        for name in re.findall(r'^\s{4}(\w+)\s*=', block, re.M):
            used = re.search(r'params\.' + name + r'\b', body) or len(re.findall(r'params\.' + name + r'\b', cfg)) > 0
            self.assertTrue(used, f'params.{name} is defined but never used')

    def test_trace_file_is_stage_specific(self):
        self.assertIn('trace.${params.stage}.txt', read('nextflow.config'))


class ResultsInternalConsistency(unittest.TestCase):
    def test_every_bench_row_sums_to_the_same_truth(self):
        for (plat, caller, regime, t, b), r in N['sv'].items():
            if b == 'ALL':
                self.assertEqual(int(r['tp_base']) + int(r['fn']), N['sv_truth_in_bed'][t], (plat, caller, regime, t))

    def test_del_plus_ins_equals_all(self):
        for plat in ('ont', 'hifi'):
            for caller in ('sniffles', 'cutesv'):
                for regime in hn.REGIMES:
                    g = lambda t, k: int(N['sv'][(plat, caller, regime, t, 'ALL')][k])
                    self.assertEqual(g('DEL', 'tp_base') + g('INS', 'tp_base'), g('ALL', 'tp_base'))

    def test_by_size_table_is_the_best_caller_default_regime(self):
        for (t, b), r in N['sv_by_size'].items():
            if t == 'ALL' and b != 'ALL':
                continue
            for p in ('ont', 'hifi'):
                src = N['sv'][(p, hn.BEST[p], 'default', t, b)]
                self.assertEqual(int(float(r[f'{p}_tp'])), int(src['tp_base']))

    def test_best_caller_is_really_best(self):
        for p in ('ont', 'hifi'):
            other = 'sniffles' if hn.BEST[p] == 'cutesv' else 'cutesv'
            self.assertGreater(hn.f1(N, p), hn.f1(N, p, caller=other))

    def test_platform_concordance_matches_truvari_summary(self):
        t = N['platform_truvari']
        self.assertEqual(int(N['platform_concordance'][('both_platforms', 'ont', 'ALL')]['n_calls']), t['TP-base'])
        self.assertEqual(int(N['platform_concordance'][('ont_only', 'ont', 'ALL')]['n_calls']), t['FN'])
        self.assertEqual(int(N['platform_concordance'][('hifi_only', 'hifi', 'ALL')]['n_calls']), t['FP'])

    def test_section12_assertions_still_hold(self):
        for p in ('ont', 'hifi'):
            self.assertGreater(float(N['qc'][p]['mapping_rate']), 0.90)
            self.assertGreater(int(N['qc'][p]['tagged_MM_reads']), 0)
            self.assertGreater(float(N['phasing'][p]['block_n50_over_read_n50']), 10)
            self.assertGreater(float(N['phasing'][p]['fraction_reads_tagged']), 0.5)
            self.assertGreater(hn.f1(N, p, field='recall'), 0.5); self.assertGreater(hn.f1(N, p, field='precision'), 0.5)
        self.assertGreater(int(N['qc']['ont']['read_n50']), int(N['qc']['hifi']['read_n50']))
        self.assertEqual(N['meth_global']['bimodal'], 'True')
        self.assertLess(float(N['meth_context'][('promoter_islands', 'promoter_cpg_island')]['mean_methylation']), 0.2)


class HeadlineNumbersInDocs(unittest.TestCase):
    """The numbers people will quote must be the numbers in results/."""

    def setUp(self):
        self.readme, self.res = squash(read('README.md')), squash(read('docs/results.md'))

    def test_imprinting_count(self):
        k = N['imprinting_pass']
        word = {3: 'three', 4: 'four', 5: 'five', 6: 'six', 7: 'seven'}[k]
        for name, doc in (('README', self.readme), ('results', self.res)):
            self.assertIn(f'{word} of seven', doc, name)
            for wrong in set(['three', 'four', 'five', 'seven']) - {word}:
                self.assertNotRegex(doc, wrong + r' of (the )?seven[^.]*(pass|cleanly)', f'{name} still claims {wrong} of seven')
        # the table verdict of every primary DMR must agree with control_pass
        table = read('docs/results.md').split('### The imprinted-locus control')[1].split('Full table')[0]
        rows = [l for l in table.splitlines() if l.startswith('| ') and not l.startswith('| DMR') and not l.startswith('|---')]
        self.assertEqual(len(rows), len(hn.PRIMARY_DMRS))
        for line, dmr in zip(rows, hn.PRIMARY_DMRS):
            verdict = line.strip().strip('|').split('|')[-1].lower()
            passed = N['imprinting'][dmr]['control_pass'] == 'True'
            self.assertEqual('yes' in verdict, passed, f'{dmr}: table says "{verdict.strip()}", control_pass={passed}')

    def test_sv_regime_table(self):
        names = {'cutesv': 'cuteSV', 'sniffles': 'Sniffles2'}
        text = read('docs/results.md').replace('**', '')
        for p, P in (('ont', 'ONT'), ('hifi', 'HiFi')):
            for c in names:
                for regime in hn.REGIMES:
                    r = N['sv'][(p, c, regime, 'ALL', 'ALL')]
                    tp, fp, fn = int(r['tp_comp']), int(r['fp']), int(r['fn'])      # from counts: the TSV stores 4 decimals and double rounding bites
                    pr, rc = tp / (tp + fp), int(r['tp_base']) / (int(r['tp_base']) + fn)
                    row = f"| {P} | {names[c]} | {regime} | {r['tp_base']} | {fp} | {fn} | {pr:.3f} | {rc:.3f} | {2*pr*rc/(pr+rc):.3f} |"
                    self.assertIn(row, text)

    def test_size_table(self):
        text = read('docs/results.md').replace('**', '')
        for t in ('DEL', 'INS'):
            for b in hn.SIZE_BINS[:4]:
                r = N['sv_by_size'][(t, b)]
                o, h = int(float(r['ont_tp'])), int(float(r['hifi_tp'])); tot = o + int(float(r['ont_fn']))
                self.assertIn(f'({o}/{tot})', text); self.assertIn(f'({h}/{tot})', text)

    def test_strict_regime_insertion_claim(self):
        d = {p: (hn.f1(N, p, svtype='INS', field='recall'), hn.f1(N, p, 'strict', svtype='INS', field='recall')) for p in ('ont', 'hifi')}
        for p in d:
            self.assertIn(f'{d[p][0]:.3f} -> {d[p][1]:.3f}', self.res)
        ont_drop, hifi_drop = (d['ont'][0] - d['ont'][1]), (d['hifi'][0] - d['hifi'][1])
        if abs(ont_drop - hifi_drop) > 0.01:
            self.assertNotIn('by the same amount', self.res)

    def test_flat_headline_values(self):
        q, ph, sm = N['qc'], N['phasing'], N['smallvar']
        must = {
            'both': [f"{int(q['ont']['read_n50']):,}", f"{int(q['hifi']['read_n50']):,}",
                     f"{hn.f1(N,'ont'):.2f}", f"{hn.f1(N,'hifi'):.2f}",
                     f"{float(sm[('ont','INDEL','*')]['f1']):.2f}", f"{float(sm[('hifi','INDEL','*')]['f1']):.2f}",
                     f"{float(q['ont']['native_coverage']):.1f}x"],
            'results': [f"{hn.f1(N,'ont'):.3f}", f"{hn.f1(N,'hifi'):.3f}", f"{N['asm_sig_sites']:,}", f"{N['asm_tested']:,}",
                        f"{len(N['asm_regions']):,} regions", f"{int(N['meth_global']['n_sites']):,}",
                        f"{int(ph['ont']['block_n50_bp'])/1e6:.2f} Mb", f"{int(ph['hifi']['block_n50_bp'])/1e6:.2f} Mb",
                        f"{100*float(ph['ont']['switch_error_rate']):.3f} %", f"{100*float(ph['hifi']['switch_error_rate']):.3f} %",
                        f"{int(ph['ont']['reads_unassigned']):,}", f"{int(ph['hifi']['reads_unassigned']):,}",
                        f"{float(q['ont']['mean_depth_target_chrom']):.1f}x", f"{float(q['hifi']['mean_depth_target_chrom']):.1f}x",
                        f"{N['sv_truth_in_bed']['ALL']} SVs", f"{N['asm_big_known']} of the top {N['asm_big_regions']}"],
        }
        for v in must['both']:
            self.assertIn(v, self.readme, f'README lacks {v}'); self.assertIn(v, self.res, f'results.md lacks {v}')
        for v in must['results']:
            self.assertIn(v, self.res, f'results.md lacks {v}')

    def test_regime_swing_ranges(self):
        rec, prec, rec_p, prec_p, missed = [], [], [], [], []
        for p in ('ont', 'hifi'):
            for c in ('cutesv', 'sniffles'):
                g = lambda regime, f: hn.f1(N, p, regime, caller=c, field=f)
                rec.append(100 * (g('default', 'recall') - g('strict', 'recall'))); prec.append(100 * (g('default', 'precision') - g('strict', 'precision')))
                rec_p.append(100 * (g('permissive', 'recall') - g('default', 'recall'))); prec_p.append(100 * (g('permissive', 'precision') - g('default', 'precision')))
                missed.append(100 * (1 - g('permissive', 'recall')))
        up = lambda x: int(x + 0.5)      # round half up (python's round() is half-to-even)
        self.assertIn(f'costs {up(min(rec))}-{up(max(rec))} points of recall and {up(min(prec))}-{up(max(prec))} points of precision', self.res)
        self.assertIn(f'{up(min(missed))}-{up(max(missed))} % of truth SVs', self.res)
        self.assertTrue(1 <= min(rec_p) and max(rec_p) <= 4.5 and 2 <= min(prec_p) and max(prec_p) <= 5.5)


class AuditArtefacts(unittest.TestCase):
    """Counts that need the uncommitted VCFs are committed as TSVs by scripts/audit_recounts.py; docs must quote them."""

    def test_truvari_containment_agrees_with_independent_recount(self):
        for t, k in (('ALL', 'all'), ('DEL', 'DEL'), ('INS', 'INS')):
            self.assertEqual(int(N['audit'][f'truth_sv_chr20_50bp_50kb_inside_bed_{k}']), N['sv_truth_in_bed'][t])

    def test_docs_quote_the_recounted_values(self):
        a = N['audit']; dec = squash(read('docs/decisions.md')); lim = squash(read('docs/limitations.md')); res = squash(read('docs/results.md'))
        self.assertIn(f"{int(a['truth_sv_chr20_50bp_50kb_all']):,}", res)
        self.assertIn(f"{a['truth_sv_chr20_50bp_50kb_DEL']} DEL, {int(a['truth_sv_chr20_50bp_50kb_INS']):,} INS", dec)
        self.assertIn(f"{int(a['truth_smvar_chr20_heterozygous']):,}", dec); self.assertIn(f"{int(a['truth_smvar_chr20_records_filter_dot']):,}", dec)
        self.assertEqual(a['truth_het_within_45kb_of_L3MBTL1_DMR'], '0')
        self.assertIn(f"{int(a['nearest_truth_het_left_of_L3MBTL1_DMR_distance_bp'])/1000:.1f} kb", lim)

    def test_pick_ac_changes_recall_not_precision_or_conclusions(self):
        for (p, c, regime), r in N['sv_ac'].items():
            base = N['sv'][(p, c, regime, 'ALL', 'ALL')]
            self.assertEqual(int(r['fp']), int(base['fp'])) if (p, c, regime) != ('hifi', 'sniffles', 'permissive') else None
            self.assertGreater(float(r['recall']), float(base['recall']) + 0.04)
        ac = lambda p, c: float(N['sv_ac'][(p, c, 'default')]['f1'])
        for p in ('ont', 'hifi'):
            other = 'sniffles' if hn.BEST[p] == 'cutesv' else 'cutesv'
            self.assertGreater(ac(p, hn.BEST[p]), ac(p, other))
        self.assertGreater(ac('ont', 'cutesv'), ac('hifi', 'sniffles'))

    def test_pick_ac_table_in_results(self):
        res = read('docs/results.md').replace('**', '')
        names = {'cutesv': 'cuteSV', 'sniffles': 'Sniffles2'}
        for p, P in (('ont', 'ONT'), ('hifi', 'HiFi')):
            for c in names:
                self.assertIn(f"| {P} | {names[c]} | {hn.f1(N, p, caller=c, field='recall'):.3f} -> {hn.ac(N, p, c):.3f} |", res)
                self.assertIn(f"-> {hn.ac(N, p, c, field='f1'):.3f} |", res)

    def test_cutesv_missed_insertion_explanation(self):
        rows = N['cutesv_missed']; res = squash(read('docs/results.md'))
        twin = sum(r['why_missed'] == 'two_truth_records_one_call' for r in rows)
        nocall = sum(r['why_missed'] == 'no_cutesv_call_within_1kb' for r in rows)
        self.assertIn(f'At {twin} of the {len(rows)}', res); self.assertIn(f'Only {nocall} of the {len(rows)}', res)
        self.assertIn(f"{sum(r['recovered_with_pick_ac'] == 'True' for r in rows)} of the {len(rows)} are recovered", res)
        self.assertNotIn('merges insertions in tandem repeats aggressively', res)


class Makefile(unittest.TestCase):
    def test_offline_targets_exist(self):
        mk = read('Makefile')
        for target in ('test:', 'numbers:'):
            self.assertIn(target, mk)


if __name__ == '__main__':
    unittest.main()
