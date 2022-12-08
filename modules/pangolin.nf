process pangolin {
    label 'pangolin'
    container = params.pangolin_docker
    conda ( (workflow.profile.contains('conda') ||  workflow.profile.contains('mamba')) ? params.pangolin_conda : null)
    publishDir "${params.output}/${params.linage_dir}/${name}", mode: params.publish_dir_mode
    
    input:
    tuple val(name), path(fasta)

    output:
    tuple val(name), path("${name}_lineage_report.csv"), emit: report

    script:
    """
    pangolin --outfile ${name}_lineage_report.csv --tempdir . --threads task.cpus ${fasta}
    """
    stub:
    """
    touch ${name}_lineage_report.csv
    """
}
