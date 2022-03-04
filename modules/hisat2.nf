process index_hisat2 {
    label 'hisat2'
    input:
    path(reference)

    output:
    tuple path(reference), path("${reference.baseName}*.ht2"), emit: index

    script:
    """
    hisat2-build -p ${task.cpus} ${reference} ${reference.baseName}
    """
    stub:
    """
    touch ${reference.baseName}.1.ht2
    """
}

process hisat2 {
    label 'hisat2'
    
    publishDir "${params.output}/${params.mapping_dir}/${name}", mode: params.publish_dir_mode

    input:
    tuple val(name), path(reads)
    tuple path(reference), path(index)

    output:
    tuple val("${name}_hisat2"), path("${name}_hisat2.bam"), emit: bam
    path "${name}_hisat2_summary.log", emit: log

    script:
    additional_parameter = params.hisat2_no_splice ? '--no-spliced-alignment' : ''
    if ( ! params.mode  == 'paired' ) {
    """
    hisat2 --rg-id=${name}_hisat2 --rg-id=PU:${name}_hisat2 --rg-id=SM:${name}_hisat2 --rg-id=PL:ILLUMINA --rg-id=LB:000 -x ${reference.baseName} -U ${reads[0]} -p ${task.cpus} --new-summary --summary-file ${name}_hisat2_summary.log ${additional_parameter} | samtools view -bS | samtools sort -o ${name}_hisat2.bam -T tmp --threads ${task.cpus}
    """
    }
    else {
    """
    hisat2 --rg-id=ID:${name}_hisat2 --rg=PU:${name}_hisat2 --rg=SM:${name}_hisat2 --rg=PL:ILLUMINA --rg=LB:000 -x ${reference.baseName} -1 ${reads[0]} -2 ${reads[1]} -p ${task.cpus} --new-summary --summary-file ${name}_hisat2_summary.log ${additional_parameter} | samtools view -bS | samtools sort -o ${name}_hisat2.bam -T tmp --threads ${task.cpus}
    """
    } 
    stub:
    """
    touch ${name}_hisat2_summary.log ${name}_hisat2.bam
    """
}