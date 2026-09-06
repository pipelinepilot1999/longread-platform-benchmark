// Clair3 (CPU) with the model matching the basecaller/chemistry, then hap.py against the v5.0q
// small-variant benchmark within its BED with homopolymer/tandem-repeat stratifications.

process CLAIR3 {
    tag "${params.platform}"
    publishDir "${params.outdir}/smallvar", mode: 'copy'

    input:
    tuple path(bam), path(bai)
    path ref
    path fai

    output:
    tuple path("${params.platform}.clair3.vcf.gz"), path("${params.platform}.clair3.vcf.gz.tbi"), emit: vcf
    path "${params.platform}.clair3.cmd"

    script:
    def model = params.clair3_model[params.platform]
    """
    export HOME=\$PWD
    cmd="/opt/bin/run_clair3.sh --bam_fn=${bam} --ref_fn=${ref} --threads=${task.cpus} --platform=${params.platform} --model_path=/opt/models/${model} --output=clair3 --ctg_name=${params.chrom} --sample_name=${params.sample}"
    echo "\$cmd" > ${params.platform}.clair3.cmd
    \$cmd
    cp clair3/merge_output.vcf.gz ${params.platform}.clair3.vcf.gz
    cp clair3/merge_output.vcf.gz.tbi ${params.platform}.clair3.vcf.gz.tbi
    """
}

process BUILD_STRATA {
    label 'analysis'
    storeDir "${params.resources}/strata_${params.chrom}"

    input:
    path scripts
    path strat_dir

    output:
    path "strata.tsv", emit: tsv
    path "*.bed.gz", emit: beds

    script:
    """
    python3 ${scripts}/build_strata.py --giab ${strat_dir} --chrom ${params.chrom} --outdir .
    """
}

process HAPPY {
    tag "${params.platform}"
    publishDir "${params.outdir}/happy", mode: 'copy'

    input:
    tuple path(vcf), path(tbi)
    tuple path(truth_vcf), path(truth_tbi), path(truth_bed)
    path ref
    path fai
    path strat_tsv
    path strat_beds

    output:
    path "${params.platform}.happy.extended.csv", emit: extended
    path "${params.platform}.happy.summary.csv", emit: summary
    path "${params.platform}.happy.*", emit: all

    script:
    """
    export HGREF=${ref}
    hap.py ${truth_vcf} ${vcf} -f ${truth_bed} -r ${ref} -l ${params.chrom} \\
        --stratification ${strat_tsv} --pass-only --threads ${task.cpus} --no-roc --no-json \\
        -o ${params.platform}.happy
    """
}
