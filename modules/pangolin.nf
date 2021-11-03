process pangolin {
    label 'pangolin'

    publishDir "${params.output}/${params.linage_dir}/${name}", mode: params.publish_dir_mode
    
    input:
    tuple val(name), path(fasta)

    output:
    tuple val(name), path("${name}_lineage_report.csv")

    script:
    """
    conda update pangolin
    pangolin --outfile ${name}_lineage_report.csv --tempdir . --threads ${task.cpus} ${fasta}
    """
}