// QC: read-length distribution (N50/median/max), read count, mapping and supplementary rates,
// mapping to chr20 versus elsewhere for the pre-selected reads, and a binned depth profile.

process READ_STATS {
    tag "${params.platform}"
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    path fastq

    output:
    path "${params.platform}.seqkit_stats.tsv", emit: stats
    path "${params.platform}.read_lengths.txt.gz", emit: lengths

    script:
    """
    seqkit stats -a -T -j ${task.cpus} ${fastq} > ${params.platform}.seqkit_stats.tsv
    seqkit fx2tab -n -l -j ${task.cpus} ${fastq} | awk '{print \$NF}' | gzip -c > ${params.platform}.read_lengths.txt.gz
    """
}

process ALIGN_QC {
    tag "${params.platform}"
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    tuple path(bam), path(bai)

    output:
    path "${params.platform}.flagstat.txt", emit: flagstat
    path "${params.platform}.idxstats.txt", emit: idxstats
    path "${params.platform}.samtools_stats.txt", emit: stats

    script:
    """
    samtools flagstat -@ ${task.cpus} ${bam} > ${params.platform}.flagstat.txt
    samtools idxstats ${bam} > ${params.platform}.idxstats.txt
    samtools stats -@ ${task.cpus} ${bam} > ${params.platform}.samtools_stats.txt
    """
}

process MOSDEPTH {
    tag "${params.platform}"
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    tuple path(bam), path(bai)

    output:
    path "${params.platform}.mosdepth.summary.txt", emit: summary
    path "${params.platform}.regions.bed.gz", emit: regions

    script:
    """
    mosdepth -t ${task.cpus} -n --fast-mode --by 10000 -c ${params.chrom} ${params.platform} ${bam}
    """
}
