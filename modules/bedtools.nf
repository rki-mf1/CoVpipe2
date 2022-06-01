// get sequence depth in single base resolution, an aggregation has to be performed at plotting level
// this data set can also be used to get the amount of bases with >0 sequence depth (genome coverage)

process get_genomecov {
    label 'bedtools'

    publishDir "${params.output}/${params.mapping_dir}/${name}", mode: params.publish_dir_mode

    input:
    tuple val(name), path(bam)

    output:
    tuple val(name), path("${name}.coverage.tsv"),        emit: tsv
    
    script:
    """
    bedtools genomecov -ibam ${bam} -d 1> ${name}.coverage.tsv
    """
    stub:
    """
    touch ${name}.coverage.tsv
    """
}

process create_low_coverage_no_del_mask {
    label 'bedtools'

    publishDir "${params.output}/${params.consensus_dir}/${name}", mode: params.publish_dir_mode

    input:
    tuple val(name), path(vcf), path(bam)

    output:
    tuple val(name), path("${name}.lowcov.bed")

    script:
    """
    bedtools genomecov -bga -ibam ${bam} | awk '\$4 < ${params.cns_min_cov}' | bedtools merge > ${name}.lowcov.bed.tmp
    bedtools subtract -a ${name}.lowcov.bed.tmp -b ${vcf} > ${name}.lowcov.bed
    """
    stub:
    """
    touch ${name}.lowcov.bed
    """
}
