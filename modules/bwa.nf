process index_bwa {
    label 'bwa'  
    
    publishDir "${params.output}", mode: 'copy', pattern: "${reference}.*"
    
    input:
    path(reference)
    
    output:
    tuple path(reference), path("${reference}.*"), emit: index

    script:
    """
    bwa index ${reference} &> /dev/null
    """
}

process bwa {
    label 'bwa'  

    input:
    tuple val(name), path(reads)
    tuple path(reference), path(index)

    output:
    tuple val(name), path("${name}.bam"), emit: bam

    script:
    """
    bwa mem -t ${task.cpus} \
        -R '@RG\\tID:${name}\\tPU:${name}\\tSM:${name}\\tPL:ILLUMINA\\  tLB:000' \
        ${reference} \
        ${reads} | \
        samtools view -Sb -@ ${task.cpus} | \
        samtools sort -@ ${task.cpus} > ${name}.bam
    """
}