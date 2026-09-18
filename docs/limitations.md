# Limitations

Final list; the results document repeats the short form.

- **One sample, one chromosome, one dataset per platform.** Platform differences are confounded with the specific library, flow cell / SMRT cell, basecaller version and run; nothing here separates "ONT versus HiFi" from "this ONT run versus this HiFi run".
- **chr20 only**, chosen for compute, and the reads were pre-selected by their original chr20 alignment in the public BAMs (ONT: Dorado/minimap2 2.24; HiFi: pbmm2 1.10 to `GRCh38_GIABv3`). Reads that only map to chr20 after re-alignment with minimap2 2.30 are therefore absent, and reads that originally mapped to chr20 but belong elsewhere are present until re-alignment moves them.
- **Methylation on ONT only** (scoping decision D-023); the Revio reads carry 5mC tags and the same analysis would run on them.
- **No basecalling** (CPU-only host); pre-basecalled reads with existing modified-base tags were used. Basecaller versions differ by platform and by date (Dorado sup v5.0.0, 2024; ccs 7.0.0 / primrose 1.4.0, 2023).
- **No genome assembly, no isoform/transcript/direct-RNA analysis, no repeat-expansion genotyping, no whole-genome run, no experimental validation, no clinical interpretation.**
- **Two SV callers out of many**; no pbsv, no SVIM, no ensemble merging.
- **Truth set is a GIAB draft (v5.0q)** and Truvari was run without `--refine` (D-010), so some complex SVs with different representations are counted as FP+FN pairs.
- **DeepVariant was not benchmarked** against Clair3 (D-011); hap.py used the xcmp engine rather than vcfeval (D-012).
- **ECR push was not performed**: the host has no AWS credentials. Images are built locally on the host from the Dockerfiles in `containers/` (see `versions.yml`).
- **Effective depth after re-alignment differs by 2x** (ONT 46.9x vs HiFi 49.0x on chr20) although both platforms entered alignment at 49.6x, because more of the pre-selected ONT reads re-align to other contigs. The direction of this bias is against ONT.
- **Half of chr20's SVs are outside the v5.0q benchmark BED** (802 of 1,620 SVs of 50 bp-50 kb are scored). Recall in the unscored, mostly repetitive half is unknown and will be lower. Only two truth deletions above 10 kb fall inside the BED, so the "ONT advantage on large events" prediction could not be tested.
- **The best caller differs by platform** (cuteSV on ONT, Sniffles2 on HiFi), so the platform comparison carries a caller difference with it; per-caller tables are given for both.
- **L3MBTL1 is unphaseable in HG002** (no heterozygous variant within 45 kb in the truth set); the imprinting control there fails for a reason the pipeline cannot fix. A trio or a second individual would resolve it.
- **NESP-AS is partially hypomethylated** in this sample (0.30 unphased), consistent with partial loss of imprinting in lymphoblastoid-line DNA; this is a sample property, not a measurement error, but it means the control is weaker there.
- **ASM candidates are unreplicated.** 8.7 % of testable CpGs are significant; the non-imprinted candidates may include haplotagging errors near block boundaries, CpG-destroying SNVs and residual 5hmC, and nothing here separates those from genuine mQTL-driven ASM.
- **ONT 6mA calls** (`A+a`) are present in the source tags and were ignored; 5hmC was not tabulated (`--modified-bases 5mC`), so a small fraction of "unmethylated" calls at 5hmC sites is expected.
- **Only one flowcell of the ONT dataset was used** (PAW70337 of two); the second would have allowed a within-platform replicate.
- **GitHub push and ECR push:** the target repository did not exist at the end of the run and the host has no AWS credentials; the repository is committed locally on the Mac and needs the user to create the remote (see the session summary).
- **The imprinting control passes at 4 of 7 loci by the rule fixed in advance**, not 5: GNAS-XL separates 0.71 vs 0.21 (0.497) and misses the 0.50 threshold. An earlier version of these docs rounded it to a pass.
- **The ONT platform Nextflow trace was overwritten by the compare stage** (both wrote the same trace file under `results/ont/pipeline_info/`); per-task ONT runtimes are not recoverable. The recorded tool commands (`*.cmd`) are intact.
- **Per-read accuracy was not measured.** No q-score or identity statistic is reported; `samtools stats` NM-based mismatch rates (ONT 4.3 %, HiFi 1.5 %; the stats files are 11-159 MB and are not committed) include true variants and SV lengths and are not an error rate.
- **Several counts cannot be re-derived from the repository alone** because the truth VCF, BAMs and call VCFs are not committed: 1,620 chr20 truth SVs of 50 bp-50 kb, 72,728 heterozygous truth records, and the zero-heterozygous-variants-within-45-kb statement for L3MBTL1. The phase-block evidence for L3MBTL1 (a 22.7 kb block built on 4 variants) is in `results/ont/phasing/ont.blocks.tsv`.
