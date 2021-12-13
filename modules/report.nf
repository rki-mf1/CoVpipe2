process desh_qc {
    label 'bcftools'

    input:
    tuple val(name), path(consensus)
    tuple val(name), path(coverage)
    val(min_cov)

    output: 
    tuple val(name), path("${name}_desh.tsv")

    script:
    """
    echo 'sample,#n,#iupac,#lowcov' > ${name}_desh.tsv
    N=\$(tail -n +2 ${consensus} | tr -cd "N" | wc -c)
    IUPAC=\$(tail -n +2 ${consensus} | tr -cd "RYSWKMBDHVN" | wc -c)
    COVTH=\$(awk -F "\t" '{{COV=\$3; if (COV < ${min_cov}) {{ print COV }} }}' ${coverage} | wc -l)
    echo ${name},\$N,\$IUPAC,\$COVTH >> ${name}_desh.tsv
    """
}

// process fastp_table {
    // label 'r'
    // label 'smallTask'

    // input:

    // output:

    // script:
    // """
    // """
// }

// process kraken_table {
    // label 'r'
    // label 'smallTask'

    // input:

    // output:

    // script:
    // """
    // """
// }

// process flagstat_table {
    // label 'r'
    // label 'smallTask'

    // input:

    // output:

    // script:
    // """
    // """
// }

// process fragment_size_table {
    // label 'r'
    // label 'smallTask'

    // input:

    // output:

    // script:
    // """
    // """
// }

// process coverage_table {
    // label 'r'
    // label 'smallTask'

    // input:

    // output:

    // script:
    // """
    // """
// }

process multiqc_report {
    label 'multiqc'
    label 'smallTask'

    input:
    path(fastp)
    path(kraken)
    path(flagstat)
    path(pangolin)

    output:
    path ("multiqc_report.html")

    script:
    """
    multiqc .
    """
}

// process rmarkdown_report {
    // label 'r'
    // label 'smallTask'

    // input:

    // output:

    // script:
    // """
    // """
// }