# Self-validation log (spec section 12)

Each assertion, what it checked, the outcome, and what was done on failure. Updated as stages complete.

## Data setup

| Assertion | Outcome |
|---|---|
| MM/ML tags present in the source BAM, verified on records | **PASS.** ONT PAW70337 `calls.sorted.bam`: 111/111 chr20 records inspected carry `MM:Z`/`ML:B` (codes `C+m`, `C+h`, and `A+a`; the 6mA calls are ignored downstream). HiFi Revio BAM: 88/100 records carry `MM:Z:C+m`/`ML` (primrose leaves reads without a confident call untagged). |
| Reference build agrees across reads, reference, both truth sets, imprinted-region annotation; contig naming consistent | **PASS.** All `chr`-prefixed GRCh38. ONT headers carry chr20 `M5:b18e6c531b0bd70e949a7fc20859cb01`; the uppercased chr20 sequence of the host FASTA has the same md5. v5.0q VCFs/BEDs, GIAB v3.6 stratifications, UCSC CpG islands, GENCODE v50 and the Akbari 2022 DMR table are GRCh38/hg38 (the DMR coordinates coincide with GENCODE gene positions for GNAS and NNAT). The Revio BAM was aligned to `GRCh38_GIABv3` (chr20 unchanged); reads are re-aligned. |
| Coverage over chr20 adequate for haplotype-partitioned analysis (flag < 20x) | **PASS.** Native chr20 coverage: ONT 54.6x, HiFi 49.6x; matched target 49.6x. |
| Both platforms are the same sample, verified from metadata | **PASS.** ONT `@RG SM:hg002 LB:hg002` (dataset `giab_2025.01/.../HG002/`); HiFi `@RG SM:HG002 LB:NIST HG002`, README names NIST RM 8392 HG002. The rejected ONT flowcell PAO83395 is the counter-example: its run metadata says `sample_id=PrePCR`. |
| Coverage-matching target computed and applied before comparison | **PASS (target).** 49.61x; ONT downsampled with fraction 0.9084 (`results/ont/qc/ont.downsample.tsv`). |

## After alignment

| Assertion | Outcome |
|---|---|
| MM/ML tags survived the FASTQ round-trip and realignment (count > 0) | **FAILED once, then PASS.** Smoke run 2026-09-05: `ASSERT_TAGS` found 0 tagged records in the realigned chr20 BAM although the FASTQ comment carried the tags. Cause: the read-group string was double-escaped in the Nextflow script (`-R '@RG\\tID:...'` reached minimap2 with literal `\\t`), and with that malformed `-R` minimap2 2.30 silently stops copying the FASTQ comment even with `-y`. Reproduced on 2 reads: double-escaped `-R` -> 0 tags, single-escaped -> MM and ML present. Fixed in `modules/minimap2_align.nf`; the assertion is what caught it, exactly as intended. Final tag-survival counts: see `results/*/qc/*.tag_survival.tsv`. |
| Mapping rate > 90 % both platforms | **PASS.** ONT 99.99 % (170,882 / 170,899 primary records), HiFi 99.94 % (191,463 / 191,582). |
| Supplementary alignments non-zero and plausible | **PASS.** ONT 16,680 supplementary records (0.098 per primary), HiFi 9,652 (0.050 per primary). The 2:1 ratio follows read length: a 30 kb read is twice as likely to cross an SV breakpoint as a 17 kb read. |
| Read N50 differs substantially, ONT longer | **PASS.** Analysed reads: ONT N50 29,418 bp (max 629 kb), HiFi N50 17,394 bp (max 42 kb). |
| Correct minimap2 preset per platform, verified from the recorded command | **PASS (from the recorded command; the ONT platform trace was later overwritten, see note below).** `results/ont/qc/ont.minimap2.cmd`: `-ax map-ont -y --secondary=no`; `results/hifi/qc/hifi.minimap2.cmd`: `-ax map-hifi -y --secondary=no`. Final tag survival: ONT 166,412 / 166,412 chr20 primary reads carry MM/ML (100 %); HiFi 166,832 / 190,001 (87.8 %, the same fraction primrose tagged in the source BAM). |

## Smoke test (50,000 ONT reads, 2026-09-05)

All 30 processes ran end to end after three fixes (missing `ps` in conda images, the minimap2 `-R` escaping above, and modkit 0.6.4 flag names). Note that `seqkit head` on a coordinate-sorted slice takes reads from the start of chr20 only (mean depth 13.9x but concentrated in the first ~14 Mb), so the smoke recall numbers (SV 0.14, SNP 0.33) measure the covered fraction, not the caller; precision (SV 0.94-0.95, SNP 0.994) and the phasing mechanics (14 blocks, N50 2.2 Mb, switch error 0.05 %) are the meaningful checks. `-resume` was exercised four times on the same session id. One more trap surfaced when moving to the full run: `--smoke false` on the command line arrives as the string `'false'`, which is truthy in Groovy, so the first "full" launch silently re-used every cached smoke task; caught because every process reported `Cached` and `results/ont` was identical to the smoke output. Fixed by comparing `params.smoke.toString() == 'true'` and by launching without the flag.

## After SV calling

| Assertion | Outcome |
|---|---|
| Both callers produced non-empty VCFs; counts within an order of magnitude | **PASS.** Inside the benchmark BED (default regime): ONT Sniffles2 594 calls, cuteSV 616; HiFi Sniffles2 590, cuteSV 500. |
| Truvari overlap between callers well above zero | **PASS.** 553 of 592 ONT Sniffles2 calls match a cuteSV call (93 %); 486 of 588 on HiFi (83 %). |
| Default-regime recall/precision not below 0.5 on either platform | **PASS.** Precision 0.92-0.97, recall 0.61-0.71 (all callers, both platforms). Recall is well above 0.5 but not high; the size-stratified tables show it is the 100-299 bp and >= 1 kb insertions that are missed, consistently on both platforms and both callers, which points at truth-set representation in tandem repeats rather than at a naming or build problem. |
| Platform-versus-platform overlap substantial | **PASS.** 539 SVs called on both platforms; 77 ONT-only, 50 HiFi-only (87.5 % of ONT cuteSV calls and 91.5 % of HiFi Sniffles2 calls in the BED are shared). |

## After phasing

| Assertion | Outcome |
|---|---|
| Phase block N50 substantially longer than read N50 | **PASS.** ONT 2.12 Mb vs 29.4 kb (72x); HiFi 0.53 Mb vs 17.4 kb (30x). |
| Haplotagged read fraction a clear majority | **PASS.** ONT 75.8 % tagged (24.2 % unassigned); HiFi 84.9 % tagged (15.2 % unassigned). See results section 3 for why ONT is the lower of the two. |

## After methylation

| Assertion | Outcome |
|---|---|
| Global CpG methylation distribution bimodal | **PASS.** 743,745 CpGs at >= 10x: 20.1 % below 0.2, 47.9 % above 0.8, 10.9 % between 0.4 and 0.6 (`results/ont/figures/fig07`). |
| Promoter CpG islands predominantly hypomethylated | **PASS.** Promoter islands: mean 0.067, 91.5 % of sites below 0.2; non-promoter islands mean 0.37; gene bodies 0.74. |
| Known imprinted DMRs intermediate unphased and cleanly separated when phased | **PARTIAL: 4 of 7 pass the pre-stated criterion; one near miss, one weak, one documented failure.** All seven primary-control DMRs are intermediate unphased (0.30-0.69). Four pass `control_pass` (MCTS2P 0.04 vs 0.99, NNAT 0.97 vs 0.29, GNAS NESP55 0.12 vs 0.97, GNAS exon 1A 0.98 vs 0.01); GNAS-XL (0.71 vs 0.21) has a separation of 0.497 and misses the 0.50 threshold, so it is a failure by the rule as written; the GNAS cluster is internally consistent: the three maternal DMRs are methylated on the same haplotype and the paternal NESP55 DMR on the other. NESP-AS separates in the right direction but weakly (0.47 vs 0.09, unphased 0.30). **L3MBTL1 does not separate (0.54 vs 0.46)**, and the diagnosis is a phasing-power failure, not a methylation failure: the DMR sits in a stretch where HG002 is essentially homozygous (the v5.0q truth and both callsets have almost no heterozygous variants within 45 kb), so the reads over the DMR are haplotagged from four weak variants in a 23 kb block and the assignment does not track parental origin. Both haplotypes therefore reproduce the unphased ~0.5, which is what a mixture looks like. On the nanopore-detected DMR windows of the same loci the verdicts are the same except NNAT, which misses the 0.3 bound there (0.303) (`results/ont/imprinting_control.tsv`). The pipeline was not modified in response; the exception is reported as a limitation of read-based phasing in a low-heterozygosity individual/region. |

## Audit note (2026-09-18)

A later audit of docs against `results/` found that the ONT trace file (kept as `results/ont/pipeline_info/trace.compare.txt`) holds only the two tasks of the compare
stage: the compare stage defaulted to `--outdir results/ont` and the trace was written with `overwrite = true`, so the ONT platform
trace was lost. The per-tool command files (`results/*/qc/*.cmd`, `sv/*.cmd`, `smallvar/*.cmd`) survive and are what the preset
assertion rests on. `nextflow.config` now writes `trace.<stage>.txt` so a compare run cannot overwrite a platform trace. The same
audit corrected the imprinting-control count from five to four of seven; `tests/test_docs_match_results.py` now enforces it.

The same audit, once the host was reachable again, re-derived every count that depends on uncommitted VCFs (all confirmed; `results/audit_truth_counts.tsv`), corrected the truth-SV type split in D-002, and tested the documented explanation for cuteSV's low HiFi insertion recall. That explanation was wrong: 22 of the 36 'missed' insertions >= 1 kb are two-record truth loci that cuteSV called once, correctly, and that `--pick single` scores as one hit and one miss. A `--pick ac` sensitivity analysis was added (D-024): recall +9 to +13 points on every caller and platform, precision and all conclusions unchanged.
