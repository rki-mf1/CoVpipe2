process index_fasta {
    label 'samtools'
    label 'smallTask'

    publishDir "${params.output}/${params.reference_dir}", mode: params.publish_dir_mode

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

process stats {
    label 'samtools'
    publishDir "${params.output}/${params.mapping_dir}/${name}", mode: params.publish_dir_mode
    
    input:
    tuple val(name), path(bam)

    output:
    tuple val(name), path("${name}.stats"), emit: stats
    tuple val(name), path("${name}.stats_small"), emit: stats_small

    script:
    """
    samtools stats ${bam} -@ $task.cpus > ${name}.stats
    grep -P 'SN\\traw total sequences:' ${name}.stats | awk -F'\\t' 'BEGIN{OFS="\\t"} {print "total", \$3}' > ${name}.stats_small
    grep -P 'SN\\treads mapped:' ${name}.stats | awk -F'\\t' 'BEGIN{OFS="\\t"} {print "mapped", \$3}' >> ${name}.stats_small
    """
    stub:
    """
    touch ${name}.stats ${name}.stats_small
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
    samtools flagstat ${bam} -@ $task.cpus > ${name}.flagstat
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
    samtools view -F 4 ${bam} -@ $task.cpus | cut -f9 1>> ${name}.fragment_size.tsv
    """
    stub:
    """
    touch ${name}.fragment_size.tsv
    """
}

process filter_isize_bam {
    label 'samtools'

    input:
    tuple val(name), path(bam)
    val(insert_size_threshold)

    output:
    tuple val(name), path("${name}_isize_filtered.bam"), emit: bam
    
    script:
    """
    samtools view -h ${bam} -@ $task.cpus | \
        awk '{ if(\$0 ~ /^@/) {print \$0} else { if(sqrt(\$9^2) < sqrt(${insert_size_threshold}^2)) {print \$0}} }' | \
        samtools view -Sbh -@ $task.cpus | \
        samtools sort -@ $task.cpus > ${name}_isize_filtered.bam
    """
    stub:
    """
    touch ${name}_isize_filtered.bam
    """
}
