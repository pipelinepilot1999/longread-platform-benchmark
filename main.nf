#!/usr/bin/env nextflow
// longread-platform-benchmark: ONT vs PacBio HiFi on HG002 chr20.
//
//   nextflow run main.nf --platform ont  --raw_bam data/ont_chr20.raw.bam  --target_cov <x>
//   nextflow run main.nf --platform hifi --raw_bam data/hifi_chr20.raw.bam --target_cov <x>
//   nextflow run main.nf --stage compare --best_caller.ont cutesv --best_caller.hifi sniffles   (after both platform runs)
//
// One platform per run is deliberate: the host cannot hold both platforms' BAMs and work directories
// at once (docs/decisions.md D-007). Arm C (methylation) runs only when --platform ont.

nextflow.enable.dsl = 2

include { BAM_TO_FASTQ; NATIVE_COVERAGE; DOWNSAMPLE; MINIMAP2_INDEX; MINIMAP2_ALIGN; SUBSET_CHROM; ASSERT_TAGS } from './modules/minimap2_align'
include { READ_STATS; ALIGN_QC; MOSDEPTH } from './modules/qc_stats'
include { SNIFFLES } from './modules/sniffles'
include { CUTESV } from './modules/cutesv'
include { TRUTH_PREP; TRUVARI_BENCH; CALLER_CONCORDANCE; PLATFORM_CONCORDANCE } from './modules/truvari_bench'
include { CLAIR3; BUILD_STRATA; HAPPY } from './modules/small_variants'
include { WHATSHAP_PHASE; WHATSHAP_STATS; WHATSHAP_COMPARE } from './modules/phasing'
include { WHATSHAP_HAPLOTAG } from './modules/haplotag'
include { MODKIT_PILEUP; MODKIT_PILEUP_PHASED } from './modules/modkit_pileup'
include { QC_SUMMARY; SV_SUMMARY; SMALLVAR_SUMMARY; PHASING_SUMMARY; METHYLATION_ANALYSIS; CROSS_PLATFORM } from './modules/asm_analysis'

workflow PLATFORM {
    if (!(params.platform in ['ont', 'hifi'])) error "--platform must be ont or hifi"
    if (!params.raw_bam) error "--raw_bam is required"

    ref     = file(params.ref_fasta, checkIfExists: true)
    fai     = file("${params.ref_fasta}.fai", checkIfExists: true)
    raw_bam = file(params.raw_bam, checkIfExists: true)
    preset  = params.mm_preset[params.platform]
    scripts = file("${projectDir}/scripts", checkIfExists: true)

    // ---- truth sets (chr20 subsets; cached in resources) ----
    TRUTH_PREP(file(params.truth_sv_vcf), file("${params.truth_sv_vcf}.tbi"), file(params.truth_sv_bed),
               file(params.truth_sm_vcf), file("${params.truth_sm_vcf}.tbi"), file(params.truth_sm_bed))

    // ---- reads -> matched coverage -> realignment ----
    BAM_TO_FASTQ(raw_bam)
    NATIVE_COVERAGE(BAM_TO_FASTQ.out.fastq)
    native_cov = NATIVE_COVERAGE.out.tsv.map { f -> f.text.readLines()[1].split('\t')[-1] as Double }
    DOWNSAMPLE(BAM_TO_FASTQ.out.fastq, native_cov)
    READ_STATS(DOWNSAMPLE.out.fastq)
    MINIMAP2_INDEX(ref, preset)
    MINIMAP2_ALIGN(DOWNSAMPLE.out.fastq, MINIMAP2_INDEX.out.mmi, preset)
    ALIGN_QC(MINIMAP2_ALIGN.out.bam)
    SUBSET_CHROM(MINIMAP2_ALIGN.out.bam)
    ASSERT_TAGS(SUBSET_CHROM.out.bam)
    MOSDEPTH(SUBSET_CHROM.out.bam)
    QC_SUMMARY(scripts, NATIVE_COVERAGE.out.tsv, DOWNSAMPLE.out.tsv, READ_STATS.out.stats, ALIGN_QC.out.flagstat,
               ALIGN_QC.out.idxstats, MOSDEPTH.out.summary, ASSERT_TAGS.out.tsv)

    // ---- Arm A: structural variants ----
    SNIFFLES(SUBSET_CHROM.out.bam, ref, fai, file(params.trf_bed))
    CUTESV(SUBSET_CHROM.out.bam, ref, fai)
    calls   = SNIFFLES.out.vcf.mix(CUTESV.out.vcf)
    regimes = Channel.fromList(params.truvari_regimes.collect { k, v -> [k, v] })
    TRUVARI_BENCH(calls.combine(regimes), TRUTH_PREP.out.sv, ref, fai)
    CALLER_CONCORDANCE(SNIFFLES.out.vcf, CUTESV.out.vcf, TRUTH_PREP.out.sv, ref, fai)
    SV_SUMMARY(scripts, TRUVARI_BENCH.out.dir.map { it[2] }.collect(), CALLER_CONCORDANCE.out.dir)

    // ---- Arm B: small variants and phasing ----
    CLAIR3(SUBSET_CHROM.out.bam, ref, fai)
    BUILD_STRATA(scripts, file(params.strat_dir))
    HAPPY(CLAIR3.out.vcf, TRUTH_PREP.out.sm, ref, fai, BUILD_STRATA.out.tsv, BUILD_STRATA.out.beds)
    SMALLVAR_SUMMARY(scripts, HAPPY.out.extended)
    WHATSHAP_PHASE(CLAIR3.out.vcf, SUBSET_CHROM.out.bam, ref, fai)
    WHATSHAP_STATS(WHATSHAP_PHASE.out.vcf)
    WHATSHAP_COMPARE(WHATSHAP_PHASE.out.vcf, TRUTH_PREP.out.sm)
    WHATSHAP_HAPLOTAG(WHATSHAP_PHASE.out.vcf, SUBSET_CHROM.out.bam, ref, fai)
    PHASING_SUMMARY(scripts, WHATSHAP_STATS.out.tsv, WHATSHAP_STATS.out.blocks, WHATSHAP_COMPARE.out.tsv,
                    WHATSHAP_HAPLOTAG.out.counts, READ_STATS.out.lengths)

    // ---- Arm C: methylation (ONT only, by design) ----
    if (params.platform == 'ont') {
        MODKIT_PILEUP(SUBSET_CHROM.out.bam, ref, fai)
        MODKIT_PILEUP_PHASED(WHATSHAP_HAPLOTAG.out.bam, ref, fai)
        METHYLATION_ANALYSIS(scripts, MODKIT_PILEUP.out.bed, MODKIT_PILEUP_PHASED.out.hp1, MODKIT_PILEUP_PHASED.out.hp2,
                             file(params.cpg_islands), file(params.gene_gtf), file(params.imprinted_bed), fai)
    }
}

// Cross-platform stage: platform concordance with Truvari, then the python comparison/cross-arm/figures script.
workflow COMPARE {
    ref = file(params.ref_fasta, checkIfExists: true)
    fai = file("${params.ref_fasta}.fai", checkIfExists: true)
    scripts = file("${projectDir}/scripts", checkIfExists: true)
    TRUTH_PREP(file(params.truth_sv_vcf), file("${params.truth_sv_vcf}.tbi"), file(params.truth_sv_bed),
               file(params.truth_sm_vcf), file("${params.truth_sm_vcf}.tbi"), file(params.truth_sm_bed))
    def best = params.best_caller ? params.best_caller : [ont: 'sniffles', hifi: 'sniffles']   // set from Arm A F1 (docs/decisions.md)
    ont  = Channel.of([best.ont,  file("results/ont/sv/ont.${best.ont}.vcf.gz"),   file("results/ont/sv/ont.${best.ont}.vcf.gz.tbi")])
    hifi = Channel.of([best.hifi, file("results/hifi/sv/hifi.${best.hifi}.vcf.gz"), file("results/hifi/sv/hifi.${best.hifi}.vcf.gz.tbi")])
    PLATFORM_CONCORDANCE(ont, hifi, TRUTH_PREP.out.sv, ref, fai)
    CROSS_PLATFORM(scripts, file('results/ont', checkIfExists: true), file('results/hifi', checkIfExists: true), PLATFORM_CONCORDANCE.out.dir, "${best.ont},${best.hifi}")
}

workflow {
    if (params.stage == 'compare') COMPARE() else PLATFORM()
}
