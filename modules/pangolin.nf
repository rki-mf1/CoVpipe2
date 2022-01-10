process pangolin {
    label 'pangolin'

    publishDir "${params.output}/${params.linage_dir}/${name}", mode: params.publish_dir_mode
    
    input:
    tuple val(name), path(fasta)
    val(pangolin_version)

    output:
    tuple val(name), path("${name}_lineage_report.csv")

    script:
    """
    pangolin --outfile ${name}_lineage_report.csv --tempdir . --threads ${task.cpus} ${fasta}
    """
}

process update_pangolin {
    // execute this loccaly - would most likely fail on custer systems
    label 'pangolin'
    executor 'local'
    cpus 1
    cache false

    output:
    env(pangolin_version)

    script:
    """
    conda update pangolin
    pangolin_version=\$(pangolin --version)
    """
}