process trim_primer {
    label 'ptrimmer' 

    publishDir "${params.output}/${params.read_dir}//${name}/primer-clipping", mode: params.publish_dir_mode

    input:
    tuple val(name), path(reads)
    path(primer_txt)

    output:
    tuple val(name), path("${name}*.noprimer.fastq.gz"),  emit: reads

    script:

    mode = params.mode == 'paired' ? 'pair' : 'single' 

    """
    ptrimmer \
        -t ${mode} \
        -a ${primer_txt} \
        -f ${reads[0]} \
        -d ${name}.R1.noprimer.fastq \
        -r ${reads[1]} \
        -e ${name}.R2.noprimer.fastq \
        -m ${params.max_primer_mismatches} \
        --keep
    gzip ${name}.R1.noprimer.fastq ${name}.R2.noprimer.fastq
    """
}
