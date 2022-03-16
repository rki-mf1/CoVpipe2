process fastqc {
    label 'fastqc'
    publishDir "${params.output}/${params.read_dir}/${name}", mode: params.publish_dir_mode

    input:
    tuple val(name), path(reads)

    output:
    tuple val(name), path("*_fastqc*"), emit: zip

    script:
    """
    fastqc -t ${task.cpus} ${reads}
    """
}