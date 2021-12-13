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