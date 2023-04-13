process nextclade {
    label 'nextclade'
    container = params.nextclade_docker
    conda ( (workflow.profile.contains('conda') ||  workflow.profile.contains('mamba')) ? params.nextclade_conda : null)
    publishDir "${params.output}/${params.linage_dir}/${name}/", mode: params.publish_dir_mode

    input:
    tuple val(name), path(consensus)
    val(nextclade_dataset_name)
    val(nextclade_dataset_version)


    output:
    tuple val(name), path("${name}_clade.tsv"), emit: results
    tuple val(name), path("*.aligned.fasta"), emit: fasta_algn
    env(used_nextclade_version), emit: version
    val(nextclade_dataset_version), emit: dataset_version

    script:
    """
    nextclade_version_curr=\$(nextclade --version)
    nextclade dataset get --name ${nextclade_dataset_name} --tag ${nextclade_dataset_version}  --output-dir 'data/${nextclade_dataset_name}_${nextclade_dataset_version}'
    nextclade run -j $task.cpus --input-dataset data/${nextclade_dataset_name}_${nextclade_dataset_version} --output-fasta ${name}.aligned.fasta --output-tsv tmp.tsv ${consensus}
    cat tmp.tsv | tr -d "\r" > ${name}_clade.tsv

    used_nextclade_version=\$nextclade_version_curr
    """
    stub:
    """
    touch ${name}_clade.tsv ${name}.aligned.fasta
    used_nextclade_version=42
    """
}
