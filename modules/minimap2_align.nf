// Read export, coverage matching and re-alignment.
//
// The chr20 slice of the public BAM is converted back to FASTQ with the modified-base tags carried
// in the comment field (samtools fastq -T MM,ML), optionally downsampled to the matched-coverage
// target, and re-aligned against the FULL reference with the platform preset and minimap2 -y so the
// comment-field tags are copied into the output BAM. Reads are then subset to chr20 for analysis.
// A naive round-trip silently destroys MM/ML (the pipeline completes and Arm C produces nothing),
// so ASSERT_TAGS counts tagged records in the aligned BAM and fails hard on zero.

process BAM_TO_FASTQ {
    tag "${params.platform}"
    label 'disk_guard'
    publishDir "${params.outdir}/qc", mode: 'copy', pattern: '*.txt'

    input:
    path raw_bam

    output:
    path "${params.platform}.raw.fastq.gz", emit: fastq
    path "${params.platform}.raw_bam.flagstat.txt", emit: flagstat

    script:
    """
    samtools flagstat -@ ${task.cpus} ${raw_bam} > ${params.platform}.raw_bam.flagstat.txt
    # -F 0x900: primary records only (supplementary/secondary carry hard-clipped partial sequence).
    # -T MM,ML,qs: keep the modified-base tags (and ONT read q-score) in the FASTQ comment field.
    samtools fastq -@ ${task.cpus} -F 0x900 -T MM,ML,qs ${raw_bam} | bgzip -@ ${task.cpus} -c > ${params.platform}.raw.fastq.gz
    n=\$(zcat ${params.platform}.raw.fastq.gz | head -400 | awk 'NR%4==1' | grep -c 'MM:Z')
    echo "reads_with_MM_in_first_100: \$n"
    """
}

process NATIVE_COVERAGE {
    tag "${params.platform}"
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    path fastq

    output:
    path "${params.platform}.native_coverage.tsv", emit: tsv

    script:
    """
    seqkit stats -a -T -j ${task.cpus} ${fastq} > seqkit.tsv
    # columns of seqkit stats -a -T: file format type num_seqs sum_len min_len avg_len max_len Q1 Q2 Q3 sum_gap N50 ...
    awk -F'\\t' -v p=${params.platform} -v L=${params.chrom_len} 'NR==1{for(i=1;i<=NF;i++)h[\$i]=i}
        NR==2{printf "platform\\treads\\tbases\\tread_n50\\tmean_len\\tnative_coverage\\n%s\\t%s\\t%s\\t%s\\t%s\\t%.3f\\n", p, \$h["num_seqs"], \$h["sum_len"], \$h["N50"], \$h["avg_len"], \$h["sum_len"]/L}' seqkit.tsv > ${params.platform}.native_coverage.tsv
    cat ${params.platform}.native_coverage.tsv
    """
}

process DOWNSAMPLE {
    tag "${params.platform}"
    label 'disk_guard'
    publishDir "${params.outdir}/qc", mode: 'copy', pattern: '*.tsv'

    input:
    path fastq
    val native_cov

    output:
    path "${params.platform}.reads.fastq.gz", emit: fastq
    path "${params.platform}.downsample.tsv", emit: tsv

    script:
    def target = params.target_cov ? params.target_cov as Double : null
    def frac   = (target && target < native_cov) ? String.format('%.5f', target / native_cov) : null
    if (params.smoke.toString() == 'true')   // a CLI '--smoke false' arrives as the string 'false', which is truthy
        """
        seqkit head -n ${params.smoke_reads} ${fastq} | bgzip -c > ${params.platform}.reads.fastq.gz
        printf 'platform\\tmode\\tnative_coverage\\ttarget_coverage\\tfraction\\n${params.platform}\\tsmoke_${params.smoke_reads}_reads\\t${native_cov}\\tNA\\tNA\\n' > ${params.platform}.downsample.tsv
        """
    else if (frac)
        """
        seqkit sample -p ${frac} -s ${params.downsample_seed} -j ${task.cpus} ${fastq} | bgzip -@ ${task.cpus} -c > ${params.platform}.reads.fastq.gz
        printf 'platform\\tmode\\tnative_coverage\\ttarget_coverage\\tfraction\\n${params.platform}\\tdownsampled\\t${native_cov}\\t${target}\\t${frac}\\n' > ${params.platform}.downsample.tsv
        """
    else
        """
        ln -s ${fastq} ${params.platform}.reads.fastq.gz
        printf 'platform\\tmode\\tnative_coverage\\ttarget_coverage\\tfraction\\n${params.platform}\\tnative\\t${native_cov}\\t${target ?: 'NA'}\\t1.0\\n' > ${params.platform}.downsample.tsv
        """
}

process MINIMAP2_INDEX {
    tag "${preset}"
    storeDir "${params.resources}/mmi"

    input:
    path ref
    val preset

    output:
    path "GRCh38.${preset}.mmi", emit: mmi

    script:
    """
    minimap2 -t ${task.cpus} -x ${preset} -d GRCh38.${preset}.mmi ${ref}
    """
}

process MINIMAP2_ALIGN {
    tag "${params.platform}"
    label 'disk_guard'
    publishDir "${params.outdir}/qc", mode: 'copy', pattern: '*.cmd'

    input:
    path fastq
    path mmi
    val preset

    output:
    tuple path("${params.platform}.full.bam"), path("${params.platform}.full.bam.bai"), emit: bam
    path "${params.platform}.minimap2.cmd", emit: cmd

    script:
    // one escaped backslash in Groovy -> literal \t reaches minimap2, which turns it into a tab
    def rg = "@RG\\tID:${params.platform}\\tSM:${params.sample}\\tPL:${params.platform == 'ont' ? 'ONT' : 'PACBIO'}"
    """
    cmd="minimap2 -t ${task.cpus} -ax ${preset} -y --secondary=no -R '${rg}' ${mmi} ${fastq}"
    echo "\$cmd" > ${params.platform}.minimap2.cmd
    minimap2 -t ${task.cpus} -ax ${preset} -y --secondary=no -R '${rg}' ${mmi} ${fastq} \\
      | samtools sort -@ 2 -m 1500M -T sort_tmp -o ${params.platform}.full.bam -
    samtools index -@ ${task.cpus} ${params.platform}.full.bam
    """
}

process SUBSET_CHROM {
    tag "${params.platform}"
    publishDir "${params.outdir}/alignment", mode: 'copy'

    input:
    tuple path(bam), path(bai)

    output:
    tuple path("${params.platform}.${params.chrom}.bam"), path("${params.platform}.${params.chrom}.bam.bai"), emit: bam

    script:
    """
    samtools view -@ ${task.cpus} -b -o ${params.platform}.${params.chrom}.bam ${bam} ${params.chrom}
    samtools index -@ ${task.cpus} ${params.platform}.${params.chrom}.bam
    """
}

process ASSERT_TAGS {
    tag "${params.platform}"
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    tuple path(bam), path(bai)

    output:
    path "${params.platform}.tag_survival.tsv", emit: tsv

    script:
    """
    total=\$(samtools view -c -F 0x904 ${bam})
    mm=\$(samtools view -c -F 0x904 -d MM ${bam})
    ml=\$(samtools view -c -F 0x904 -d ML ${bam})
    printf 'platform\\tprimary_reads\\twith_MM\\twith_ML\\tfraction_MM\\n' > ${params.platform}.tag_survival.tsv
    awk -v p=${params.platform} -v t=\$total -v a=\$mm -v b=\$ml 'BEGIN{printf "%s\\t%d\\t%d\\t%d\\t%.4f\\n", p, t, a, b, (t>0?a/t:0)}' >> ${params.platform}.tag_survival.tsv
    cat ${params.platform}.tag_survival.tsv
    if [ "\$mm" -eq 0 ] || [ "\$ml" -eq 0 ]; then
        echo "ASSERTION FAILED: no MM/ML-tagged records survived realignment (section 4.3 handling broke)" >&2
        exit 1
    fi
    """
}
