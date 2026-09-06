// Assign each read to haplotype 1, 2 or unassigned with the phased VCF; the unassigned fraction is a
// direct measure of phasing power. Supplementary alignments inherit their primary's tag.
process WHATSHAP_HAPLOTAG {
    tag "${params.platform}"
    publishDir "${params.outdir}/phasing", mode: 'copy', pattern: '*.{tsv,tsv.gz}'
    publishDir "${params.outdir}/alignment", mode: 'copy', pattern: '*.bam*'

    input:
    tuple path(vcf), path(tbi)
    tuple path(bam), path(bai)
    path ref
    path fai

    output:
    tuple path("${params.platform}.${params.chrom}.haplotagged.bam"), path("${params.platform}.${params.chrom}.haplotagged.bam.bai"), emit: bam
    path "${params.platform}.haplotag_list.tsv.gz", emit: list
    path "${params.platform}.haplotag_counts.tsv", emit: counts

    script:
    """
    whatshap haplotag --reference ${ref} --ignore-read-groups --tag-supplementary \\
        --output-threads ${task.cpus} --output-haplotag-list ${params.platform}.haplotag_list.tsv.gz \\
        -o ${params.platform}.${params.chrom}.haplotagged.bam ${vcf} ${bam}
    samtools index -@ ${task.cpus} ${params.platform}.${params.chrom}.haplotagged.bam
    total=\$(samtools view -c -F 0x904 ${params.platform}.${params.chrom}.haplotagged.bam)
    h1=\$(samtools view -c -F 0x904 -d HP:1 ${params.platform}.${params.chrom}.haplotagged.bam)
    h2=\$(samtools view -c -F 0x904 -d HP:2 ${params.platform}.${params.chrom}.haplotagged.bam)
    printf 'platform\\tprimary_reads\\thp1\\thp2\\tunassigned\\tfraction_tagged\\n' > ${params.platform}.haplotag_counts.tsv
    awk -v p=${params.platform} -v t=\$total -v a=\$h1 -v b=\$h2 'BEGIN{printf "%s\\t%d\\t%d\\t%d\\t%d\\t%.4f\\n", p, t, a, b, t-a-b, (t>0?(a+b)/t:0)}' >> ${params.platform}.haplotag_counts.tsv
    cat ${params.platform}.haplotag_counts.tsv
    """
}
