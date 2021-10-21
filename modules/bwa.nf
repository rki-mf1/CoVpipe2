process index_bwa {
    label 'bwa'  
    
    publishDir "${params.output}", mode: 'copy', pattern: "${name}_bwa_index.tar.gz"
    
    input:
    tuple val(name), path(reference)
    
    output:
    tuple path(reference), path("${name}_bwa_index.tar.gz"), emit: index

    script:
    """
    bwa index ${reference} &> /dev/null
    tar -zcvf ${name}_bwa_index.tar.gz ${reference}.*
    """
}

process bwa {
    label 'bwa'  

    input:
    tuple val(name), path(reads)
    tuple path(reference), path(index)

    script:
    """
    (   time \
                bwa mem -t ${task.cpus} \
                    -R '@RG\tID:${name}\tPU:${name}\tSM:${name}\tPL:ILLUMINA\tLB:000' \
                    ${reference} \
                    ${reads} | \
                    samtools view -Sb -@ ${task.cpus} -o ${name}.bam \
            ) &> ${name}_map2reference.log
    """
}