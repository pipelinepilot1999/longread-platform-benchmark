# longread-platform-benchmark

**Oxford Nanopore vs PacBio HiFi on the same sample (HG002), the same region (chr20), through one pipeline:
structural variants, read-based phasing, and haplotype-resolved methylation with an imprinting positive control.**

The deliverable is the evaluation, not the pipeline: every number in `docs/results.md` comes with the
matching regime, coverage and truth set that produced it, and every parameter choice is logged in
`docs/decisions.md`.

**Summary.** On HG002 chr20 at a matched 49.6x, Oxford Nanopore (R10.4.1, Dorado sup v5.0.0) and PacBio Revio HiFi
reach almost the same structural-variant accuracy against GIAB v5.0q (best-caller F1 0.81 vs 0.79, Truvari default regime),
and the expected split by SV class only half appears: ONT recall is equal or higher in every insertion size bin and in every
deletion bin above 100 bp, HiFi wins only 50-99 bp deletions (0.77 vs 0.71) and precision (0.94 vs 0.93). Where the platforms
really separate is small indels and phasing: HiFi indel F1 0.87 vs ONT 0.75 (in homopolymers >= 12 bp, 0.71 vs 0.40) while
SNV F1 is 0.996 vs 0.993, and ONT read length turns into 2.1 Mb phase blocks against 0.53 Mb for HiFi with a lower switch
error (0.07 % vs 0.22 %). ONT phased methylation recovers the imprinting signature at five of seven known chr20 DMRs cleanly
(GNAS cluster internally consistent: maternal DMRs on one haplotype, the paternal NESP55 DMR on the other); the one clear
failure, L3MBTL1, is a run of homozygosity where read-based phasing has nothing to phase.

## Scope, stated plainly

- **chr20 only.** This is a compute constraint (8 vCPU, CPU-only, shared host), not a methods claim. Nothing here is a whole-genome analysis.
- **Reads were pre-selected by their original chr20 alignment** in the public BAMs (ONT: Dorado's minimap2 2.24; HiFi: pbmm2 1.10), converted back to FASTQ and re-aligned against the full GRCh38 with minimap2 2.30. Reads that would only reach chr20 after re-alignment are excluded; reads that originally landed on chr20 but belong elsewhere are removed by the re-alignment (the fraction that stays on chr20 is reported).
- **Same individual on both platforms:** HG002 (NIST RM 8391/8392, Coriell NA24385). Sample identity was checked from read-group / run metadata, not filenames.
- **Coverage matching:** native chr20 coverage was computed for each platform (total bases of primary chr20 reads / 64,444,167 bp); the deeper platform was downsampled to the shallower one on the FASTQ before alignment (`seqkit sample`, fixed seed). Native and matched coverages are in `results/qc_summary.tsv`.
- **Methylation was run on ONT only, by design.** The Revio reads carry 5mC tags too; the ONT modified-base ecosystem is the more established one and running the arm twice adds disk and complexity without a new finding (D-014). It is a scoping decision, not an oversight.
- **No basecalling.** Pre-basecalled reads with existing modified-base tags were used on both platforms.
- **One dataset per platform**, so platform differences are confounded with the specific library, flow cell / SMRT cell, basecaller and run.

## Data

| | ONT | PacBio HiFi |
|---|---|---|
| Source | `ont-open-data` `giab_2025.01`, HG002 flowcell PAW70337 | GIAB `PacBio_HiFi-Revio_20231031`, HG002 48x GRCh38-GIABv3 BAM |
| Instrument / chemistry | PromethION 48, FLO-PRO114M (R10.4.1), SQK-LSK114, run 2024-06-19 | Revio, SMRTbell prep kit 3.0, two SMRT cells (2023-09/10), 24 h movies |
| Basecaller / modified bases | Dorado `dna_r10.4.1_e8.2_400bps_sup@v5.0.0`, `5mC_5hmC@v2.0.1` | ccs 7.0.0, primrose 1.4.0 (5mC) |
| Sample identity (metadata) | `@RG SM:hg002 LB:hg002` | `@RG SM:HG002 LB:NIST HG002`, README: RM 8392 |
| Native chr20 coverage | 54.6x | 49.6x |
| Matched coverage | 49.6x | 49.6x |
| Read N50 (analysed reads) | {{ONT_N50}} bp | {{HIFI_N50}} bp |

Reference: `GCA_000001405.15_GRCh38_no_alt_analysis_set` (the chr20 MD5 matches the ONT alignment headers).
Truth: GIAB **v5.0q** for SVs (`stvar`) and for phased small variants (`smvar`), both GRCh38, both from the
T2T-HG002-Q100 v1.1 assembly; GIAB v3.6 stratifications. v5.0q is a GIAB *draft* release. Imprinted-DMR
annotation: Akbari et al. 2022 (eLife 11:e77898), hg38.

## Pipeline

Nextflow DSL2, one Docker image per tool, every version pinned (`versions.yml`, `containers/`).
One platform per run, then a comparison entry:

```bash
scripts/run_platform.sh ont  data/ont_chr20.raw.bam  --target_cov 49.6
scripts/run_platform.sh hifi data/hifi_chr20.raw.bam --target_cov 49.6
nextflow run main.nf --stage compare --best_caller.ont cutesv --best_caller.hifi sniffles
```

```
main.nf, nextflow.config
modules/   minimap2_align  qc_stats  sniffles  cutesv  truvari_bench  small_variants  phasing  haplotag  modkit_pileup  asm_analysis
containers/ Dockerfile.{minimap2,sniffles,cutesv,truvari,variantcall,happy,phasing,modkit,analysis}
scripts/   build_strata.py qc_summary.py sv_summary.py smallvar_summary.py phasing_summary.py methylation_analysis.py platform_compare.py
results/   qc_summary.tsv sv_benchmark.tsv sv_concordance.tsv sv_platform_by_size.tsv smallvar_benchmark.tsv phasing_stats.tsv
           methylation_by_context.tsv asm_sites.tsv imprinting_control.tsv cross_arm.tsv figures/
docs/      results.md decisions.md limitations.md
```

### The MM/ML tag-preservation handling

A naive BAM -> FASTQ -> BAM round-trip silently drops the modified-base tags: `samtools fastq` discards
them and minimap2 ignores unknown tags. The pipeline uses `samtools fastq -T MM,ML,qs` to carry the tags
in the FASTQ comment and `minimap2 -y` to copy them into the output BAM, then **asserts** that tagged
records exist in the aligned chr20 BAM (`ASSERT_TAGS`, fails the run on zero). Tag survival:
ONT 166,412 / 166,412 chr20 primary reads (100 %); HiFi 166,832 / 190,001 (87.8 %, identical to the fraction primrose tagged in the source BAM).

### Arm A: structural variants
Sniffles2 2.8.0 and cuteSV 2.1.4 on the same BAM (documented per-platform settings, min length 50 bp),
Truvari 5.4.0 against v5.0q inside the benchmark BED under **three regimes**:

| Regime | `--refdist` | `--pctseq` | `--pctsize` | type | measures |
|---|---|---|---|---|---|
| permissive | 1000 | 0 | 0.5 | ignored | upper bound on recall |
| default | 500 | 0.7 | 0.7 | must match | comparability with published work (Truvari defaults) |
| strict | 100 | 0.9 | 0.9 | must match | conservative lower bound |

### Arm B: small variants and phasing
Clair3 2.0.2 (model matched to the basecaller: `r1041_e82_400bps_sup_v500` / `hifi_revio`), hap.py vs
v5.0q with homopolymer stratifications, WhatsHap 2.8 phasing/haplotagging, switch error vs the assembly-phased truth.

### Arm C: methylation (ONT)
modkit 0.6.4 CpG pileups, unphased and partitioned by the WhatsHap haplotag; global distribution, genomic
context, per-CpG allele-specific methylation (Fisher exact, BH-FDR), and the **imprinted-locus control**:
GNAS cluster, L3MBTL1, NNAT and MCTS2P must be ~50 % methylated unphased and split cleanly by haplotype.
Result: five of seven primary-control DMRs pass cleanly, NESP-AS separates weakly in the right direction, and L3MBTL1 cannot be separated because HG002 carries no heterozygous variant within 45 kb of it (see `docs/results.md` section 4 and `docs/validation.md`).

## What this project does not do
No basecalling or signal-level work; no methylation on HiFi (ONT-only by design); no genome assembly; no
isoform, transcript or direct-RNA analysis; no repeat-expansion genotyping; no whole-genome analysis; no
clinical or diagnostic interpretation; no experimental validation of any call. Two SV callers of many.
Images were built on the EC2 host and not pushed to ECR (no credentials on the host).

## Reproducing
Environment: EC2 r5.2xlarge (8 vCPU, 62 GB), Docker, Nextflow 26.04. Slice the public BAMs
(`samtools view -b -X <bam_url> <bai_url> chr20`), place the reference and resources as in
`nextflow.config`, build images with `scripts/build_containers.sh`, run the three commands above.
