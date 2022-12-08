process index_bwa {
    label 'bwa'
    
    input:
    path(reference)
    
    output:
    tuple path(reference), path("${reference}.*"), emit: index

    script:
    """
    echo $task.cpus
    bwa index ${reference} &> /dev/null
    """
    stub:
    """
    touch ${reference}.amb ${reference}.ann ${reference}.bwt ${reference}.pac ${reference}.sa
    """
}

process bwa {
    label 'bwa'

    publishDir "${params.output}/${params.mapping_dir}/${name}", mode: params.publish_dir_mode

    input:
    tuple val(name), path(reads)
    tuple path(reference), path(index)

    output:
    tuple val(name), path("${name}.bam"), emit: bam

    script:
    """
    bwa mem -t $task.cpus \
        -R '@RG\\tID:${name}\\tPU:${name}\\tSM:${name}\\tPL:ILLUMINA\\tLB:000' \
        ${reference} \
        ${reads} | \
        samtools view -Sb -@ $task.cpus | \
        samtools sort -@ $task.cpus > ${name}.bam
    """
    stub:
    """
    touch ${name}.bam
    """
}