process president {
    label 'president'

    publishDir "${params.output}/${params.report_dir}/${name}/", pattern: "${name}*.fasta", mode: params.publish_dir_mode

    input:
    tuple val(name), path(fasta)
    path(reference_fasta)
    
    output:
    tuple val(name), path("${name}_valid.fasta"), emit: valid
    tuple val(name), path("${name}_invalid.fasta"), emit: invalid
    tuple val(name), path("${name}_report.tsv"), emit: report
    
    script:
    """
    president -r ${reference_fasta} -t ${task.cpus} -q ${fasta} -x 0.90 -n 0.05 -p . -f ${name}_
    """
    stub:
    """
    touch ${name}_valid.fasta ${name}_invalid.fasta ${name}_report.tsv
    """
}