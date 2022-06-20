process pangolin {
    label 'pangolin'
    container = params.pangolin_docker
    if ( workflow.profile.contains('conda') ||  workflow.profile.contains('mamba') ) { conda = params.pangolin_conda }
    publishDir "${params.output}/${params.linage_dir}/${name}", mode: params.publish_dir_mode
    
    input:
    tuple val(name), path(fasta)

    output:
    tuple val(name), path("${name}_lineage_report.csv"), emit: report

    script:
    def args = params.pangolin_scorpio ? '' : '--skip-scorpio'
    """
    pangolin ${args} --outfile ${name}_lineage_report.csv --tempdir . --threads ${task.cpus} ${fasta}
    """
    stub:
    """
    touch ${name}_lineage_report.csv
    used_pangolin_version=42
    scorpio_version=42
    scorpio_constellations_version=42
    """
}
