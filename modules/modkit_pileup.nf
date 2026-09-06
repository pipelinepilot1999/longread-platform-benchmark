// modkit: per-read MM/ML -> per-CpG methylation frequency. 5mC only (--modified-bases 5mC; the 5hmC and
// 6mA calls in the ONT tags are not tabulated), both strands of a CpG combined. The phased pileup
// (--phased) splits records on the HP tag written by haplotagging into hp1/hp2/combined tables.
process MODKIT_PILEUP {
    tag "${params.platform}"
    publishDir "${params.outdir}/methylation", mode: 'copy'

    input:
    tuple path(bam), path(bai)
    path ref
    path fai

    output:
    path "${params.platform}.cpg.unphased.bed.gz", emit: bed
    path "${params.platform}.modkit.unphased.log"

    script:
    """
    modkit pileup ${bam} unphased.bed --ref ${ref} --cpg --combine-strands --modified-bases 5mC \\
        --region ${params.chrom} --threads ${task.cpus} --log-filepath ${params.platform}.modkit.unphased.log
    bgzip -@ ${task.cpus} -c unphased.bed > ${params.platform}.cpg.unphased.bed.gz
    """
}

process MODKIT_PILEUP_PHASED {
    tag "${params.platform}"
    publishDir "${params.outdir}/methylation", mode: 'copy'

    input:
    tuple path(bam), path(bai)
    path ref
    path fai

    output:
    path "${params.platform}.cpg.hp1.bed.gz", emit: hp1
    path "${params.platform}.cpg.hp2.bed.gz", emit: hp2
    path "${params.platform}.cpg.combined.bed.gz", emit: combined
    path "${params.platform}.modkit.phased.log"

    script:
    """
    mkdir -p phased
    modkit pileup ${bam} phased --phased --prefix cpg --ref ${ref} --cpg --combine-strands --modified-bases 5mC \\
        --region ${params.chrom} --threads ${task.cpus} --log-filepath ${params.platform}.modkit.phased.log
    ls -l phased
    bgzip -@ ${task.cpus} -c \$(ls phased/*hp1*) > ${params.platform}.cpg.hp1.bed.gz
    bgzip -@ ${task.cpus} -c \$(ls phased/*hp2*) > ${params.platform}.cpg.hp2.bed.gz
    bgzip -@ ${task.cpus} -c \$(ls phased/*combined*) > ${params.platform}.cpg.combined.bed.gz
    """
}
