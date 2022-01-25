process pangolin {
    label 'pangolin'

    publishDir "${params.output}/${params.linage_dir}/${name}", mode: params.publish_dir_mode
    
    input:
    tuple val(name), path(fasta)

    output:
    tuple val(name), path("${name}_lineage_report.csv"), emit: report
    env(pangolin_version), emit: version

    script:
    """
    pangolin --outfile ${name}_lineage_report.csv --tempdir . --threads ${task.cpus} ${fasta}
    pangolin_version=\$(pangolin --version)
    """
}

process update_pangolin {
    // execute this loccaly - would most likely fail on custer systems
    label 'pangolin'
    executor 'local'
    cpus 1
    cache false

    script:
    """
    conda update pangolin
    """
}