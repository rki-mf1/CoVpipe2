process index_fasta {
    label 'samtools'  

    input:
    path(reference)

    output:
    path("${reference}.fai"), emit: index

    script:
    """
    samtools faidx ${reference}
    """
}

process index_bam {
    label 'samtools'  

    input:
    tuple val(name), path(bam)

    output:
    tuple val(name), path("${bam}.bai"), emit: index

    script:
    """
    samtools index ${bam}
    """
}
