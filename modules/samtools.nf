process index_samtools {

    label 'samtools'  

    publishDir "${params.output}", mode: 'copy', pattern: "*.fai"

    input:
    tuple val(name), file(reference)

    output:
    tuple val(name), file("${reference}.fai"), emit: index

    script:
    """
    samtools faidx ${reference}
    """
}
