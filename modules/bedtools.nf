// get sequence depth in single base resolution, an aggregation has to be performed at plotting level
// this data set can also be used to get the amount of bases with >0 sequence depth (genome coverage)

process get_genomecov {
    label 'bedtools'

    input:
    tuple val(name), path(bam)

    output:
    tuple val(name), path("${name}.coverage.tsv"),        emit: tsv
    tuple val(name), path("${name}_getCoverage.log"),     emit: log
    
    script:
    """
    (bedtools genomecov -ibam ${bam} -d 1> ${name}.coverage.tsv) 2> ${name}_getCoverage.log
    """
}
