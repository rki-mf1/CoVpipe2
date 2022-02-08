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
    stub:
    """
    touch ${reference}.fai
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
    stub:
    """
    touch ${bam}.bai
    """
}

process flagstat {
    label 'samtools'
    
    input:
    tuple val(name), path(bam)

    output:
    tuple val(name), path("${name}.flagstat"), emit: flagstat
    tuple val(name), path("${name}.flagstat.csv"), emit: csv

    script:
    """
    samtools flagstat ${bam} -@ ${task.cpus} > ${name}.flagstat
    cat ${name}.flagstat | sed -e 's/ + /;/' | sed -e 's/ /;/' 1> ${name}.flagstat.csv
    """
    stub:
    """
    touch ${name}.flagstat ${name}.flagstat.csv
    """
}

process get_fragment_size {
    label 'samtools'
    
    input:
    tuple val(name), path(bam)

    output:
    tuple val(name), path("${name}.fragment_size.tsv")

    script:
    """
    echo 0 1> ${name}_fragment_size.tsv
    samtools view -F 4 ${bam} -@ ${task.cpus} | cut -f9 1>> ${name}.fragment_size.tsv
    """
    stub:
    """
    touch ${name}.fragment_size.tsv
    """
}