// Read-based phasing of Clair3 PASS calls with WhatsHap (no reference panel), block statistics,
// and switch/Hamming error against the assembly-phased v5.0q truth.

process WHATSHAP_PHASE {
    tag "${params.platform}"
    publishDir "${params.outdir}/phasing", mode: 'copy'

    input:
    tuple path(vcf), path(tbi)
    tuple path(bam), path(bai)
    path ref
    path fai

    output:
    tuple path("${params.platform}.phased.vcf.gz"), path("${params.platform}.phased.vcf.gz.tbi"), emit: vcf

    script:
    """
    bcftools view -f PASS -r ${params.chrom} -Oz -o pass.vcf.gz ${vcf}
    tabix -p vcf pass.vcf.gz
    whatshap phase --reference ${ref} --ignore-read-groups -o ${params.platform}.phased.vcf.gz pass.vcf.gz ${bam}
    tabix -p vcf ${params.platform}.phased.vcf.gz
    """
}

process WHATSHAP_STATS {
    tag "${params.platform}"
    publishDir "${params.outdir}/phasing", mode: 'copy'

    input:
    tuple path(vcf), path(tbi)

    output:
    path "${params.platform}.whatshap_stats.tsv", emit: tsv
    path "${params.platform}.blocks.tsv", emit: blocks
    path "${params.platform}.blocks.gtf", emit: gtf

    script:
    """
    whatshap stats --tsv ${params.platform}.whatshap_stats.tsv --block-list ${params.platform}.blocks.tsv --gtf ${params.platform}.blocks.gtf ${vcf}
    """
}

process WHATSHAP_COMPARE {
    tag "${params.platform}"
    publishDir "${params.outdir}/phasing", mode: 'copy'

    input:
    tuple path(vcf), path(tbi)
    tuple path(truth_vcf), path(truth_tbi), path(truth_bed)

    output:
    path "${params.platform}.whatshap_compare.tsv", emit: tsv
    path "${params.platform}.switch_errors.bed", emit: bed

    script:
    """
    whatshap compare --names truth,${params.platform} --ignore-sample-name --only-snvs \\
        --tsv-pairwise ${params.platform}.whatshap_compare.tsv --switch-error-bed ${params.platform}.switch_errors.bed \\
        ${truth_vcf} ${vcf} > ${params.platform}.whatshap_compare.log
    """
}
