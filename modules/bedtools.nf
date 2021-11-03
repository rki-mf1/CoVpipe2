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
}