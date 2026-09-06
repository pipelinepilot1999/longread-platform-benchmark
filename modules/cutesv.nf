// cuteSV with the platform-specific clustering parameters from its README (its defaults are not
// appropriate across platforms); genotyping on; min SV length 50 bp.
process CUTESV {
    tag "${params.platform}"
    publishDir "${params.outdir}/sv", mode: 'copy'

    input:
    tuple path(bam), path(bai)
    path ref
    path fai

    output:
    tuple val('cutesv'), path("${params.platform}.cutesv.vcf.gz"), path("${params.platform}.cutesv.vcf.gz.tbi"), emit: vcf
    path "${params.platform}.cutesv.cmd"

    script:
    def args = params.cutesv_args[params.platform]
    """
    mkdir -p work
    cmd="cuteSV -t ${task.cpus} -l ${params.min_sv_len} --genotype -S ${params.sample}_${params.platform} ${args} ${bam} ${ref} raw.vcf work"
    echo "\$cmd" > ${params.platform}.cutesv.cmd
    \$cmd
    bcftools sort -Oz -o ${params.platform}.cutesv.vcf.gz raw.vcf
    tabix -p vcf ${params.platform}.cutesv.vcf.gz
    """
}
