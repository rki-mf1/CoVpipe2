process index_picard {

    label 'gatk'  

    publishDir "${params.output}", mode: 'copy', pattern: "*.dict"
    publishDir "${params.output}", mode: 'copy', pattern: "*.log"

    input:
    tuple val(name), file(reference)

    output:
    tuple val(name), file("*.dict"), emit: index
    tuple val(name), file("*.log"),  emit: log

    script:
    """
    gatk CreateSequenceDictionary -R ${reference} &> ${name}.indexPicard.log
    """
}

