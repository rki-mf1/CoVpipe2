process index_bowtie2 {
    label 'bowtie2'

    input: 
    path(reference)

    output:
    tuple path(reference), path("${reference.baseName}*.bt2")

    """
    bowtie2-build -f ${reference} ${reference.baseName}
    """
}

process bowtie2 {
    label 'bowtie2'

    input:
    tuple val(name), path(reads)
    tuple path(reference), path(index)

    output:
    tuple val("${name}_bowtie2"), path("${name}_bowtie2.bam"), emit: bam

    script:
    """
    bowtie2 --rg-id=${name}_bowtie2 --rg=PU:${name}_bowtie2 --rg=SM:${name}_bowtie2 --rg=PL:ILLUMINA --rg=LB:000 --end-to-end -p ${task.cpus} --very-sensitive -X 1000 -x ${reference.baseName} -1 ${reads[0]} -2 ${reads[1]} | samtools view -Sb -@ ${task.cpus} | \
        samtools sort -@ ${task.cpus} > ${name}_bowtie2.bam
    """
}