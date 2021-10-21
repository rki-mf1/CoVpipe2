process dos2unix {
    label 'dos2unix'  

    publishDir "${params.output}", mode: 'copy', pattern: "${name}.prepared.fasta"

    input:
    path(reference)

    output:
    path("${reference.simpleName}.prepared.fasta"), emit: reference

    script:
    """
    dos2unix ${reference}
    mv ${reference} ${reference.simpleName}.prepared.fasta
    """
}
