process dos2unix {
    label 'dos2unix'  

    publishDir "${params.output}", mode: 'copy', pattern: "${name}.prepared.fasta"

    input:
    tuple val(name), path(reference)

    output:
    tuple val(name), path("${name}.prepared.fasta"), emit: reference

    script:
    """
    dos2unix ${reference}
    mv ${reference} ${name}.prepared.fasta
    """
}
