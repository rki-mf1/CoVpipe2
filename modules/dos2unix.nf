process dos2unix {
    label 'dos2unix'

    publishDir "${params.output}/${params.reference_dir}", mode: params.publish_dir_mode

    input:
    path(reference)

    output:
    path("${reference.simpleName}.prepared.fasta"), emit: reference

    script:
    """
    dos2unix ${reference}
    mv ${reference} ${reference.simpleName}.prepared.fasta
    """
    stub:
    """
    touch ${reference.simpleName}.prepared.fasta
    """
}
