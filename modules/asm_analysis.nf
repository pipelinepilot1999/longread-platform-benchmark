// Python analysis steps (container: analysis). Per-platform summaries plus the methylation / ASM /
// imprinting-control analysis (ONT only), and the cross-platform + cross-arm comparison stage.

process QC_SUMMARY {
    label 'analysis'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path scripts
    path native_cov
    path downsample
    path seqkit
    path flagstat
    path idxstats
    path mosdepth_summary
    path tag_survival

    output:
    path "qc_summary.${params.platform}.tsv", emit: tsv

    script:
    """
    python3 ${scripts}/qc_summary.py --platform ${params.platform} --native ${native_cov} --downsample ${downsample} \\
        --seqkit ${seqkit} --flagstat ${flagstat} --idxstats ${idxstats} --mosdepth ${mosdepth_summary} --tags ${tag_survival} \\
        --chrom ${params.chrom} --out qc_summary.${params.platform}.tsv
    """
}

process SV_SUMMARY {
    label 'analysis'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path scripts
    path bench_dirs
    path concordance_dir

    output:
    path "sv_benchmark.${params.platform}.tsv", emit: bench
    path "sv_concordance.${params.platform}.tsv", emit: concordance
    path "sv_calls_annotated.${params.platform}.tsv.gz", emit: calls

    script:
    """
    python3 ${scripts}/sv_summary.py --platform ${params.platform} --bench ${bench_dirs} --concordance ${concordance_dir} --outdir .
    """
}

process SMALLVAR_SUMMARY {
    label 'analysis'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path scripts
    path extended

    output:
    path "smallvar_benchmark.${params.platform}.tsv", emit: tsv

    script:
    """
    python3 ${scripts}/smallvar_summary.py --platform ${params.platform} --extended ${extended} --out smallvar_benchmark.${params.platform}.tsv
    """
}

process PHASING_SUMMARY {
    label 'analysis'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path scripts
    path stats
    path blocks
    path compare
    path haplotag_counts
    path read_lengths

    output:
    path "phasing_stats.${params.platform}.tsv", emit: tsv

    script:
    """
    python3 ${scripts}/phasing_summary.py --platform ${params.platform} --stats ${stats} --blocks ${blocks} \\
        --compare ${compare} --haplotag ${haplotag_counts} --read-lengths ${read_lengths} --out phasing_stats.${params.platform}.tsv
    """
}

process METHYLATION_ANALYSIS {
    label 'analysis'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path scripts
    path unphased
    path hp1
    path hp2
    path cpg_islands
    path gtf
    path imprinted_bed
    path fai

    output:
    path "methylation_by_context.tsv", emit: context
    path "asm_sites.tsv", emit: asm
    path "asm_regions.tsv", emit: asm_regions
    path "imprinting_control.tsv", emit: imprinting
    path "methylation_global.tsv", emit: global
    path "methylation_site_counts.tsv", emit: counts
    path "figures/*.png", emit: figures

    script:
    """
    mkdir -p figures
    python3 ${scripts}/methylation_analysis.py --unphased ${unphased} --hp1 ${hp1} --hp2 ${hp2} \\
        --cpg-islands ${cpg_islands} --gtf ${gtf} --imprinted ${imprinted_bed} --fai ${fai} --chrom ${params.chrom} \\
        --min-cov ${params.min_site_cov_unphased} --min-cov-hap ${params.min_site_cov} --outdir .
    """
}

process CROSS_PLATFORM {
    label 'analysis'
    publishDir "results", mode: 'copy'

    input:
    path scripts
    path ont_dir
    path hifi_dir
    path concordance_dir
    val best

    output:
    path "*.tsv", emit: tables
    path "figures/*.png", emit: figures

    script:
    """
    mkdir -p figures
    python3 ${scripts}/platform_compare.py --ont ${ont_dir} --hifi ${hifi_dir} --concordance ${concordance_dir} \\
        --chrom ${params.chrom} --chrom-len ${params.chrom_len} --best ${best} --outdir .
    """
}
