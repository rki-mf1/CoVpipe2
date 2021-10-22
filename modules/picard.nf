process index_picard {
    label 'gatk'

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

