process index_fasta {
    label 'samtools'
    label 'smallTask'

    publishDir "${params.output}/${params.genome_dir}", mode: params.publish_dir_mode

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
    label 'smallTask'
    
    publishDir "${params.output}/${params.mapping_dir}/${name}", mode: params.publish_dir_mode

    input:
    tuple val(name), path(bam)

    output:
    tuple val(name), path (bam), path("${bam}.bai")

    script:
    """
    samtools index ${bam}
    """
}
