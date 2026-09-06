// Truth-set preparation and Truvari benchmarking.
//
// TRUTH_PREP subsets both v5.0q benchmarks to chr20 and removes ALT="*" records from the SV truth
// (v5.0q README: truvari bench mis-categorises them). TRUVARI_BENCH runs one caller under one
// matching regime, restricted to the benchmark BED. CALLER_CONCORDANCE compares the two callers
// against each other inside the same BED so every shared/unique call can be cross-referenced to
// its truth status from the default-regime bench. PLATFORM_CONCORDANCE does the same across platforms.

process TRUTH_PREP {
    storeDir "${params.resources}/truth_${params.chrom}"

    input:
    path sv_vcf
    path sv_tbi
    path sv_bed
    path sm_vcf
    path sm_tbi
    path sm_bed

    output:
    tuple path("truth_sv.${params.chrom}.vcf.gz"), path("truth_sv.${params.chrom}.vcf.gz.tbi"), path("truth_sv.${params.chrom}.bed"), emit: sv
    tuple path("truth_sm.${params.chrom}.vcf.gz"), path("truth_sm.${params.chrom}.vcf.gz.tbi"), path("truth_sm.${params.chrom}.bed"), emit: sm

    script:
    """
    bcftools view -r ${params.chrom} -e 'ALT="*"' -Oz -o truth_sv.${params.chrom}.vcf.gz ${sv_vcf}
    tabix -p vcf truth_sv.${params.chrom}.vcf.gz
    awk '\$1=="${params.chrom}"' ${sv_bed} > truth_sv.${params.chrom}.bed
    # v5.0q flags records where the Q100 assembly is itself heterozygous/gapped (HET1/HET2/GAP1/GAP2); keep only
    # unflagged records so hap.py and whatshap compare are scored against unambiguous truth (D-003).
    bcftools view -r ${params.chrom} -f .,PASS -Oz -o truth_sm.${params.chrom}.vcf.gz ${sm_vcf}
    tabix -p vcf truth_sm.${params.chrom}.vcf.gz
    awk '\$1=="${params.chrom}"' ${sm_bed} > truth_sm.${params.chrom}.bed
    """
}

process TRUVARI_BENCH {
    tag "${params.platform}:${caller}:${regime}"
    publishDir "${params.outdir}/truvari", mode: 'copy'

    input:
    tuple val(caller), path(vcf), path(tbi), val(regime), val(args)
    tuple path(truth_vcf), path(truth_tbi), path(truth_bed)
    path ref
    path fai

    output:
    tuple val(caller), val(regime), path("${caller}_${regime}"), emit: dir

    script:
    """
    truvari bench -b ${truth_vcf} -c ${vcf} -f ${ref} --includebed ${truth_bed} \\
        --passonly --sizemin ${params.min_sv_len} --sizefilt 30 ${args} -o ${caller}_${regime}
    echo "truvari bench --includebed ${truth_bed} --passonly --sizemin ${params.min_sv_len} --sizefilt 30 ${args}" > ${caller}_${regime}/cmd.txt
    """
}

process CALLER_CONCORDANCE {
    tag "${params.platform}"
    publishDir "${params.outdir}/truvari", mode: 'copy'

    input:
    tuple val(c1), path(v1), path(t1)
    tuple val(c2), path(v2), path(t2)
    tuple path(truth_vcf), path(truth_tbi), path(truth_bed)
    path ref
    path fai

    output:
    path "caller_concordance", emit: dir

    script:
    """
    truvari bench -b ${v1} -c ${v2} -f ${ref} --includebed ${truth_bed} \\
        --passonly --sizemin ${params.min_sv_len} --sizefilt 30 ${params.truvari_regimes['default']} -o caller_concordance
    echo "base=${c1} comp=${c2} regime=default" > caller_concordance/cmd.txt
    """
}

process PLATFORM_CONCORDANCE {
    tag "${c1}_ont_vs_${c2}_hifi"
    publishDir "${params.outdir}/truvari", mode: 'copy'

    input:
    tuple val(c1), path(v1), path(t1)
    tuple val(c2), path(v2), path(t2)
    tuple path(truth_vcf), path(truth_tbi), path(truth_bed)
    path ref
    path fai

    output:
    path "platform_concordance", emit: dir

    script:
    """
    truvari bench -b ${v1} -c ${v2} -f ${ref} --includebed ${truth_bed} \\
        --passonly --sizemin ${params.min_sv_len} --sizefilt 30 ${params.truvari_regimes['default']} -o platform_concordance
    echo "base=ont:${c1} comp=hifi:${c2} regime=default" > platform_concordance/cmd.txt
    """
}
