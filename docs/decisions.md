# Decision log

One entry per decision where a reasonable alternative existed. Numbered in the order the decision was
taken. "Effect on results" was filled in once the relevant stage had run; no entry is
updated at the end of the run.

### D-001: Reference build and FASTA
**Chose:** GRCh38, `GCA_000001405.15_GRCh38_no_alt_analysis_set` (195 contigs, `chr`-prefixed), the copy already on the host at `/home/ec2-user/wgs-wes/ref/GRCh38.fa` (md5 `a6da8681616c05eb542f1d91606a7b2f`).
**Alternatives considered:** GRCh37 (the only build for the older NIST SV v0.6 benchmark); CHM13v2.0 (v5.0q is released for it, and the Revio BAM exists for it, but no ONT alignment or GIAB stratifications are available in the same form); the GIAB `GRCh38_GIABv3` masked reference used for the Revio BAM.
**Reasoning:** every input agrees on this build: the ONT BAM/CRAM headers carry `M5:b18e6c531b0bd70e949a7fc20859cb01` for chr20, identical to the md5 of the uppercased chr20 sequence in the host FASTA; the v5.0q benchmarks, GIAB v3.6 stratifications, UCSC CpG islands and GENCODE v50 are all GRCh38. The Revio BAM was aligned to `GRCh38_GIABv3` (no-alt set with masked false duplications on other chromosomes plus decoys); chr20 is unchanged in that reference, and the reads are re-aligned here anyway, so only the chr20 pre-selection depends on it.
**Source:** file headers inspected with `samtools view -H`; GIAB v5.0q README.
**Effect on results:** none directly; it is what makes the SV, small-variant and imprinting annotations comparable.

### D-002: Truth sets
**Chose:** GIAB **v5.0q** (2025-01-17, T2T-HG002-Q100 v1.1 assembly-based) for both arms: `HG002_GRCh38_v5.0q_stvar` (SV) and `HG002_GRCh38_v5.0q_smvar` (small variants, fully phased).
**Alternatives considered:** `NIST_SV_v0.6` (Tier 1, 2019) — GRCh37 only, so it cannot be used with GRCh38 reads without lifting over. `CMRG_v1.00` SV benchmark — GRCh38 but restricted to 273 medically relevant genes (few on chr20). `NISTv4.2.1` small variants — GRCh38, mapping-based, the copy on the host carries no phasing (0 `|` genotypes on chr20), so it cannot score switch errors. Q100 `T2TQ100-V1.1` — not present on the GIAB FTP any more; v5.0q is its successor.
**Reasoning:** v5.0q is the only current GIAB release that provides SV and small-variant benchmarks on the same build from the same underlying assembly, and its small-variant VCF is phased chromosome-wide (assembly haplotypes), which is exactly what switch-error scoring needs. It is labelled a *draft* by GIAB; that is recorded in the limitations.
**Source:** GIAB FTP listing surveyed 2026-09-04; `NIST_HG002_v5.0q_variant-benchmarksets_README.md`.
**Effect on results:** chr20 SV truth after filtering: 610 DEL and 1,102 INS of >= 50 bp (mostly INS 300-999 bp); benchmark BED covers 58.6 Mb of chr20 (214 regions).

### D-003: Truth-set filtering
**Chose:** SV truth: chr20 records with `ALT="*"` removed (11,595 on chr20). Small-variant truth: chr20 records with FILTER `.` only; records flagged `HET1/HET2/GAP1/GAP2` (assembly heterozygous/gapped, ~60k of 180k chr20 records) dropped.
**Alternatives considered:** use the VCFs as released and rely on the benchmark BED alone; keep flagged records.
**Reasoning:** the v5.0q README states truvari bench mis-categorises `ALT=*` records and recommends removing them. The HET/GAP flags mark places where the Q100 assembly itself is ambiguous; keeping them would score calls against uncertain truth and would inject unphased/inconsistent genotypes into the switch-error comparison. The BED already excludes most of them.
**Source:** v5.0q README; own judgment for the FILTER rule.
**Effect on results:** chr20 small-variant truth after the FILTER rule: 120,027 records (72,728 heterozygous); 60k flagged records dropped. SV truth inside the BED: 802 SVs of 50 bp-50 kb (254 DEL, 548 INS), half of the 1,620 such SVs on chr20; the excluded half is not scored and this is stated in the results. The with/without comparison was not run (time).

### D-004: ONT dataset
**Chose:** ONT open data `giab_2025.01`, HG002 flowcell **PAW70337**, `basecalling/sup/HG002/PAW70337/calls.sorted.bam`: PromethION 48 (PC48LILITH), FLO-PRO114M R10.4.1, SQK-LSK114, run 2024-06-19, Dorado `dna_r10.4.1_e8.2_400bps_sup@v5.0.0` with modified-base model `5mC_5hmC@v2.0.1`. Read group `SM:hg002 LB:hg002`. MM/ML present on 111/111 inspected chr20 records. Yield 152.8 Gb basecalled (about 49x raw).
**Alternatives considered:** (a) `giab_2023.05` PAO89685 — same sample (run metadata `sample_id=NA_24385`), R10.4.1/LSK114 but basecalled with Guppy 6.5.5 (older model, materially lower raw accuracy) and the CRAM has no `@RG`, so sample identity rests on the run summary only. (b) `giab_2023.05` PAO83395 — its own run metadata labels the sample `PrePCR`/`PreCR_restart`, not NA24385; rejected on the "verify from metadata" rule. (c) GIAB `UCSC_Ultralong_OxfordNanopore_Promethion` (2020) — R9.4.1, Guppy 3.x, no modified-base tags. (d) second 2025.01 flowcell PAW71238 — same library and model; one flowcell already gives ~49x, so merging would only add download and alignment time.
**Reasoning:** most recent model, exact matching Clair3 model available (`r1041_e82_400bps_sup_v500`), tags verified on records, sample recorded in the read group, and native depth close to the HiFi set so coverage matching costs little.
**Source:** S3 listings, `final_summary`/`throughput` files, BAM headers and records inspected on 2026-09-04.
**Effect on results:** native chr20 coverage recorded in `results/qc_summary.tsv`.

### D-005: PacBio HiFi dataset
**Chose:** GIAB `PacBio_HiFi-Revio_20231031/HG002_PacBio-HiFi-Revio_20231031_48x_GRCh38-GIABv3.bam`: Revio, two SMRT cells (m84039_230928_213653_s3, m84039_231005_222902_s1), SMRTbell prep kit 3.0, ccs 7.0.0, 5mC calls by primrose 1.4.0 (`MM:Z:C+m`), `@RG SM:HG002`, 48x, N50 17.2 kb.
**Alternatives considered:** `PacBio_CCS_15kb_20kb_chemistry2` (Sequel II, 2019, ccs 4.0.0, ~52x merged) — older chemistry/ccs, no methylation tags, and the merged BAM is a haplotagged 10x-phased file; `pacbcloud` Revio public datasets — not reachable from the host at the time.
**Reasoning:** current instrument and chemistry, same NIST reference material (RM 8392) as the ONT run, sample name in the read group, indexed over HTTPS for remote slicing.
**Source:** `README_HG002-PacBio-Revio.md`; BAM header and records.
**Effect on results:** native chr20 coverage recorded in `results/qc_summary.tsv`.

### D-006: Coverage matching
**Chose:** native coverage of each platform = total bases of primary chr20-preselected reads / 64,444,167; target = the smaller of the two; the deeper platform is downsampled to that target with `seqkit sample -p <target/native> -s 11` on the FASTQ before alignment.
**Alternatives considered:** downsampling the BAM after alignment (would leave the alignment step at unmatched depth); matching on mosdepth mean depth after alignment (circular: the downsampling has to happen before the alignment it is measured on); no matching (measures depth, not platform).
**Reasoning:** both platforms end up at the same number of sequenced bases over chr20 before any caller runs, which is what "matched coverage" means for SV recall.
**Source:** own judgment.
**Effect on results:** ONT 54.6x native, downsampled by 0.9084 to 49.6x; HiFi 49.6x native (fraction 0.99996). After re-alignment to the full genome the effective chr20 depths are 46.9x (ONT) and 49.0x (HiFi) because 9 % of the pre-selected ONT reads (3 % of HiFi) re-align elsewhere; this 2x gap works against ONT in every recall comparison and is reported alongside the matched target.

### D-007: One platform per run; host budget
**Chose:** the Nextflow pipeline processes one platform per invocation (`--platform ont|hifi`) and the comparison is a separate stage (`--stage compare`, D-022). Executor envelope 8 cpus / 30 GB; alignment and Clair3 capped at 6 cpus; a disk guard aborts any BAM-writing task below 40 GB free.
**Alternatives considered:** a single run over both platforms with `maxForks`; full use of the 8 vCPU.
**Reasoning:** the host is shared with two other running pipelines (load ~15 on 8 vCPU at start; 116 GB free of 500 GB with other projects growing), so sequential platforms with intermediate cleanup and a conservative envelope are the only way to keep all three moving. Spec §2 asks for exactly this sequencing.
**Source:** `uptime`, `df` and `docker ps` on the host, 2026-09-04.
**Effect on results:** none on numbers; wall-clock only.

### D-008: Re-alignment
**Chose:** `minimap2 2.30`, `-ax map-ont` / `-ax map-hifi`, `-y` (copy FASTQ comment tags MM/ML/qs into the BAM), `--secondary=no`, read group per platform, prebuilt full-genome `.mmi` per preset; sorted with samtools 1.21; then subset to chr20.
**Alternatives considered:** using the public alignments directly (ONT: minimap2 2.24 inside Dorado; HiFi: pbmm2 1.10) — would make the platform comparison partly an aligner comparison and the HiFi BAM carries HiPhase haplotags that must not leak into Arm B; aligning to a chr20-only FASTA — forces reads from elsewhere to mismap onto chr20 (spec §4.3); keeping secondary alignments.
**Reasoning:** one aligner, one version, platform presets, tags preserved; secondaries are not used by either SV caller or Clair3 and only inflate file size.
**Source:** minimap2 manual; spec §4.3.
**Effect on results:** mapping and supplementary rates and the fraction of pre-selected reads landing back on chr20 are in `results/qc_summary.tsv`.

### D-009: SV caller parameters
**Chose:** Sniffles2 2.8.0: defaults, `--minsvlen 50`, `--tandem-repeats human_GRCh38_no_alt_analysis_set.trf.bed` (the annotation shipped in the Sniffles repository; improves repeat-mediated calls), same settings on both platforms (Sniffles2 has no per-platform preset). cuteSV 2.1.4: `-l 50 --genotype` plus the README-recommended clustering parameters — ONT `--max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3`; HiFi `--max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5`.
**Alternatives considered:** cuteSV defaults (tuned for neither platform); Sniffles `--mosaic`; a minimum SV length of 30 bp.
**Reasoning:** the documented recommendations are the only defensible choice for a benchmark; 50 bp is the GIAB/Truvari convention for "structural".
**Source:** Sniffles2 and cuteSV READMEs.
**Effect on results:** calls inside the BED (default regime): ONT Sniffles2 594 / cuteSV 616, HiFi Sniffles2 590 / cuteSV 500. The cuteSV HiFi parameters cost recall (0.61 vs 0.69 for Sniffles2 on the same BAM), mostly insertions >= 1 kb in tandem repeats; the ONT parameters made cuteSV the better ONT caller (F1 0.807 vs 0.781).

### D-010: Truvari regimes
**Chose:** three regimes, all with `--includebed <v5.0q stvar BED> --passonly --sizemin 50 --sizefilt 30` and the default `--sizemax 50000`:
permissive `--refdist 1000 --pctseq 0 --pctsize 0.5 --typeignore`; default `--refdist 500 --pctseq 0.7 --pctsize 0.7` (Truvari 5 defaults, written explicitly); strict `--refdist 100 --pctseq 0.9 --pctsize 0.9`. No `--refine`.
**Alternatives considered:** the v5.0q README's `--refine -r 2000 -C 5000 --pick ac` recipe; `--pick multi`; dropping `--pctseq` for symbolic deletions.
**Reasoning:** the three regimes bracket the answer (upper bound on recall, comparability with published work, conservative floor). `--refine` re-represents variants through alignment, which makes per-call TP/FP status untraceable to the original VCFs, and the caller/platform concordance analyses depend on that traceability; it is noted as a follow-up in the limitations. Truvari uses the reference to compare symbolic `<DEL>` sequences, so `--pctseq` applies to both callers' deletions.
**Source:** Truvari 5 documentation; v5.0q README.
**Effect on results:** permissive -> default -> strict moves F1 0.839 -> 0.807 -> 0.733 (ONT cuteSV) and 0.828 -> 0.793 -> 0.738 (HiFi Sniffles2); the ONT-HiFi difference is +1.1 / +1.4 / -0.5 F1 points, so the platform conclusion is regime-independent. Without `--refine`, insertions in tandem repeats with shifted breakpoints or < 70-90 % allele identity are scored as FP+FN pairs; that is the main reason strict costs 7-14 precision points.

### D-011: Small-variant caller
**Chose:** Clair3 2.0.2 (CPU), model `r1041_e82_400bps_sup_v500` for ONT (exact match to the basecaller model) and `hifi_revio` for HiFi; PASS calls only downstream.
**Alternatives considered:** DeepVariant (ONT_R104 / PACBIO models) on CPU. Not run: a CPU DeepVariant call of chr20 at ~45x is several hours per platform on the 6 cores available here, on a host already loaded by two other pipelines, and DeepVariant's ONT model is trained on R10.4.1 Dorado v4/v5 sup data without an exact model per basecaller, whereas Clair3 ships one for this basecaller. The brief benchmark the spec asks for was therefore replaced by this reasoning and is listed in the limitations.
**Source:** Clair3 model list in the container (`/opt/models`); own judgment on runtime.
**Effect on results:** SNV F1 0.9934 (ONT) / 0.9957 (HiFi); indel F1 0.747 / 0.867; indel F1 inside homopolymers >= 12 bp 0.395 / 0.712. Clair3 ran in 1.5-2 h per platform on 6 cpus.

### D-012: Small-variant benchmarking
**Chose:** hap.py 0.3.15 with its built-in `xcmp` engine, `--pass-only`, restricted to chr20 (`-l chr20`) and the v5.0q smvar BED, with GIAB v3.6 stratifications subset to chr20: homopolymers 4-6 / 7-11 / >=12 bp, imperfect homopolymers >=11, all vs not-in homopolymers >=7, all vs not-in tandem repeats + homopolymers, segmental duplications.
**Alternatives considered:** the `vcfeval` engine recommended by the v5.0q README (rtg-tools is not in the hap.py biocontainer); vcfdist.
**Reasoning:** xcmp is the hap.py default and the same engine is used for both platforms, so the platform comparison is unaffected; absolute numbers may differ slightly from a vcfeval run.
**Source:** v5.0q README; hap.py documentation.
**Effect on results:** absolute numbers may differ slightly from a vcfeval run; both platforms were scored with the same engine and settings so the comparison is unaffected. The stratification is what exposes the homopolymer gap (D-011 numbers).

### D-013: Phasing
**Chose:** WhatsHap 2.8: `whatshap phase --reference --ignore-read-groups` on Clair3 PASS calls (SNVs and indels), `whatshap stats` for blocks, `whatshap compare --only-snvs --ignore-sample-name` against the v5.0q smvar truth (assembly-phased, so a single chromosome-wide haplotype pair) for switch and Hamming error, `whatshap haplotag --tag-supplementary` for read assignment.
**Alternatives considered:** LongPhase 2.0.2 (faster, SV-aware); Clair3's built-in `--enable_phasing`; HiPhase (HiFi only).
**Reasoning:** WhatsHap provides the block statistics and the compare tool in one place and behaves identically on both platforms; LongPhase would have added a second phasing model to explain. Switch error on SNVs only keeps indel representation differences out of the comparison.
**Source:** WhatsHap documentation.
**Effect on results:** ONT 109 blocks, N50 2.12 Mb, 91.8 % of hets phased, switch error 0.071 %, 75.8 % of reads tagged; HiFi 194 blocks, N50 0.53 Mb, 95.7 % phased, 0.217 %, 84.9 % tagged. `--only-snvs` keeps indel representation differences out of the switch-error count.

### D-014: Methylation pileup and coverage thresholds
**Chose:** modkit 0.6.4 `pileup --cpg --combine-strands --modified-bases 5mC` (5mC only, both strands of each CpG combined; the 5hmC and 6mA calls present in the Dorado tags are not tabulated), default pass-probability filtering; haplotype partition with `--phased` on the `HP` tag written by WhatsHap. Minimum coverage 10 reads per CpG for unphased analyses, 5 reads per haplotype for the phased comparison.
**Alternatives considered:** keeping 5hmC as a separate modification; per-strand sites; thresholds of 3 or 10 per haplotype. (`--ignore h`, the flag used in older modkit releases, no longer exists in 0.6.4; the smoke test caught this.)
**Reasoning:** the analysis is about 5mC at CpGs; combining strands doubles per-site depth. At ~45x total, after haplotagging each haplotype carries roughly 15-20 reads at a typical CpG, so 5 per haplotype retains most sites while excluding the tail where a Fisher test cannot reach significance anyway. `results/methylation_site_counts.tsv` reports sites surviving 1/3/5/10/15/20.
**Source:** modkit documentation; own judgment.
**Effect on results:** CpGs testable on both haplotypes: 732,441 at >= 1 read, 710,713 at >= 5 (chosen), 640,464 at >= 10, 446,225 at >= 15, 161,663 at >= 20. Unphased sites at >= 10: 743,745.

### D-015: Allele-specific methylation test
**Chose:** per CpG with >=5 reads on each haplotype: Fisher exact test on (modified, canonical) counts hp1 vs hp2; Benjamini-Hochberg across all tested sites on chr20; "significant" = q < 0.05 and |hp1 - hp2| >= 0.5. Regions = runs of significant sites with the same sign, gap <= 500 bp, >= 5 CpGs. Regions overlapping a known imprinted DMR are labelled controls; everything else is a *candidate*.
**Alternatives considered:** DSS/methylKit-style region tests (need replicates); a two-proportion z-test; effect-size threshold 0.3.
**Reasoning:** Fisher is exact at low counts; the 0.5 effect size targets the all-or-none pattern of imprinting and mono-allelic methylation rather than small shifts that are more likely haplotagging noise.
**Source:** own judgment.
**Effect on results:** 61,671 significant CpGs (8.7 % of tested), 1,222 regions of >= 5 CpGs; 7 of the 42 regions with >= 20 significant CpGs are known imprinted DMRs. With a 0.3 effect size the count would be larger and less specific; the candidates are reported as candidates only.

### D-016: Imprinted-DMR annotation and the control criterion
**Chose:** known imprinted DMRs on chr20 from Akbari et al. 2022 (eLife 11:e77898), Supplementary file 1 (DMRs reported by Court 2014, Joshi 2016, Zink 2018, Hernandez Mora 2018; hg38), plus their Supplementary file 2 (nanopore-detected DMRs) as a secondary annotation. Primary controls: MCTS2P/HM13, NNAT/BLCAP, L3MBTL1, and the GNAS cluster (NESP55 [paternal], NESP-AS/GNAS-AS1, GNAS-XL, GNAS exon 1A [maternal]). A DMR passes if unphased methylation is 0.25-0.75 over >= 5 CpGs and, phased, the haplotypes are separated by >= 0.5 with one haplotype >= 0.7 and the other <= 0.3.
**Alternatives considered:** deriving windows from imprinted-gene promoters (self-consistent but circular); Court et al. 2014 hg19 coordinates lifted over.
**Reasoning:** the eLife table is already on hg38 (its GNAS and NNAT coordinates coincide with the GENCODE v50 gene positions on the host), cites the primary literature per DMR, and gives the parental origin, which lets the report check that the maternal DMRs of the GNAS cluster fall on the same haplotype as each other and opposite the paternal NESP55 DMR.
**Source:** https://elifesciences.org/articles/77898 (Supplementary files 1 and 2).
**Effect on results:** 5 of 7 primary DMRs pass, NESP-AS weak (right direction, 0.47 vs 0.09), L3MBTL1 fails because HG002 has no heterozygous variant within 45 kb of it (zero in v5.0q); the criterion was not relaxed to make it pass and the failure is reported with its cause.

### D-017: Genomic-context definitions
**Chose:** CpG islands = UCSC `cpgIslandExt` (hg38); shores = 2 kb flanks of islands not in an island; promoters = TSS +/- 1 kb of GENCODE v50 basic protein-coding genes; gene bodies = protein-coding gene spans; intergenic = everything else. Sites are assigned exclusively in that priority order. Promoter islands are islands overlapping a promoter window.
**Alternatives considered:** Ensembl regulatory build; +/-2.5 kb promoters; all gene types.
**Reasoning:** conventional definitions with a known expected answer (promoter islands hypomethylated), computed from annotations already on the same build.
**Source:** own judgment.
**Effect on results:** promoter CpG islands mean 0.067 (91.5 % of sites < 0.2), non-promoter islands 0.37, gene bodies 0.74, intergenic 0.61; island methylation rises monotonically with distance to TSS (0.08 within 1 kb to 0.58 beyond 100 kb).

### D-018: Tool versions and containers
**Chose:** one Dockerfile per tool, built on the EC2 host from pinned conda packages (micromamba 1.5.10 base) or pinned upstream images: minimap2 2.30, samtools 1.21, seqkit 2.13.0, Sniffles 2.8.0, cuteSV 2.1.4, Truvari 5.4.0, Clair3 2.0.2 (`hkubal/clair3:v2.0.2`), hap.py 0.3.15 (`quay.io/biocontainers/hap.py:0.3.15--py27hcb73b3d_0`), WhatsHap 2.8, modkit 0.6.4, mosdepth 0.3.14 (`quay.io/biocontainers/mosdepth:0.3.14--h87be163_2`), bcftools/htslib 1.21, analysis image python 3.12 with pinned numpy/pandas/scipy/statsmodels/matplotlib/pysam.
**Alternatives considered:** single fat image; bare biocontainer images (most lack bcftools/samtools/tabix needed in the same step).
**Reasoning:** reproducibility with a short per-tool environment; ECR push needs AWS credentials the host does not have (see limitations).
**Source:** quay.io tag listings on 2026-09-04.
**Effect on results:** none.

### D-019: Nextflow caching and session handling
**Chose:** `process.cache = 'deep'` from the first run; runs launched with `setsid` from a wrapper script and resumed by explicit session id (`-resume <uuid>`); `MINIMAP2_INDEX`, `TRUTH_PREP` and `BUILD_STRATA` use `storeDir` under the shared resources directory so both platform runs reuse them.
**Alternatives considered:** default lenient/timestamp cache; `-resume` without a session id.
**Reasoning:** the default cache hashes name/size/timestamp and silently invalidates everything below a regenerated identical file; `-resume` picks the most recent session, which after any preview or smoke run is the wrong one. Both failure modes are silent.
**Source:** prior experience on this host.
**Effect on results:** none.

### D-020: Running the two platform passes concurrently
**Chose:** after the ONT full run started, the HiFi run was launched alongside it (separate Nextflow session, same work directory, 8-cpu envelope each).
**Alternatives considered:** strictly sequential passes as in D-007 and spec section 2.
**Reasoning:** the sequential rule exists for disk (200 GB in the spec). By the time the full runs started the volume had been grown to 750 GB with 329 GB free, and the two neighbouring pipelines had finished, so the host was idle. Two concurrent runs oversubscribe the 8 vCPU during alignment and Clair3 but overlap the many I/O-bound and single-threaded steps, which shortens total wall-clock. Peak memory stays below the 62 GB of RAM (16 GB + 14 GB per run at most), with 64 GB of swap behind it.
**Source:** `df`, `uptime`, `free` on the host, 2026-09-05.
**Effect on results:** none on numbers (each run is its own session and output directory); wall-clock only.

### D-021: Best caller per platform for the cross-platform comparison
**Chose:** cuteSV for ONT (default-regime F1 0.807 vs 0.781) and Sniffles2 for HiFi (0.793 vs 0.748), selected on the default regime inside the benchmark BED.
**Alternatives considered:** the same caller on both platforms (Sniffles2 on both would have given ONT 0.781 vs HiFi 0.793 and reversed the headline sign in the ALL row); the union or intersection of both callers.
**Reasoning:** spec section 6 asks for "the better-performing caller on each platform". Using one caller for both would fold a caller effect into the platform comparison; the per-caller tables are reported for both so a reader can make the other choice. The union/intersection would be a third pipeline, not a platform comparison.
**Source:** `results/sv_benchmark.tsv`.
**Effect on results:** the size-stratified comparison and the platform concordance use ONT cuteSV vs HiFi Sniffles2; with Sniffles2 on both, the deletion crossover and the insertion pattern are unchanged in direction but the ONT margins shrink by 2-5 points.

### D-022: Named stages instead of `-entry`
**Chose:** `--stage platform|compare` dispatched from the entry workflow.
**Alternatives considered:** `-entry COMPARE` as written in the spec's execution plan.
**Reasoning:** Nextflow 26.04 with the strict parser rejects `-entry` ("use a param to run a named workflow from the entry workflow"); the failure was immediate and the fix is a one-line dispatch.
**Source:** Nextflow error message on the first compare launch, 2026-09-05.
**Effect on results:** none.
