// Sniffles2 with its documented settings: default clustering for either platform, the bundled GRCh38
// tandem-repeat annotation (improves repeat-mediated calls), and min SV length 50 bp.
process SNIFFLES {
    tag "${params.platform}"
    publishDir "${params.outdir}/sv", mode: 'copy'

    input:
    tuple path(bam), path(bai)
    path ref
    path fai
    path trf_bed

    output:
    tuple val('sniffles'), path("${params.platform}.sniffles.vcf.gz"), path("${params.platform}.sniffles.vcf.gz.tbi"), emit: vcf
    path "${params.platform}.sniffles.cmd"

    script:
    """
    cmd="sniffles --input ${bam} --vcf raw.vcf --reference ${ref} --tandem-repeats ${trf_bed} --minsvlen ${params.min_sv_len} --threads ${task.cpus} --sample-id ${params.sample}_${params.platform}"
    echo "\$cmd" > ${params.platform}.sniffles.cmd
    \$cmd
    bcftools sort -Oz -o ${params.platform}.sniffles.vcf.gz raw.vcf
    tabix -p vcf ${params.platform}.sniffles.vcf.gz
    """
}
