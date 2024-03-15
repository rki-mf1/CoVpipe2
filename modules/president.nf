process president {
    label 'president'

    publishDir "${params.output}/${params.consensus_dir}/${name}/", pattern: "${name}*_report.tsv", mode: params.publish_dir_mode

    input:
    tuple val(name), path(fasta)
    path(reference_fasta)
    val(seq_threshold)
    val(n_threshold)
    
    output:
    tuple val(name), path("${name}_report.tsv"), path("${name}_valid.fasta"), emit: valid
    tuple val(name), path("${name}_report.tsv"), path("${name}_invalid.fasta"), emit: invalid
        
    script:
    """
    president -r ${reference_fasta} -t $task.cpus -q ${fasta} -x ${seq_threshold} -n ${n_threshold} -p . -f ${name}_
    """
    stub:
    """
    touch ${name}_valid.fasta ${name}_invalid.fasta ${name}_report.tsv
    """
}