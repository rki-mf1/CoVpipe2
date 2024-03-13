process nextclade {
    label 'nextclade'
    conda params.nextclade_conda
    container params.nextclade_docker
    publishDir "${params.output}/${params.linage_dir}/${name}/", mode: params.publish_dir_mode

    input:
    tuple val(name), path(consensus)
    val(nextclade_dataset_name)
    val(nextclade_dataset_tag)


    output:
    tuple val(name), path("${name}_clade.tsv"), emit: results
    tuple val(name), path("*.aligned.fasta"), emit: fasta_algn
    env(used_nextclade_version), emit: version
    env(used_nextcladedataset_info), emit: dataset_info

    script:
    def tag = params.update ? '' : "--tag ${nextclade_dataset_tag}"
    """
    nextclade_version_curr=\$(nextclade --version)
    nextclade dataset get --name ${nextclade_dataset_name} ${tag} --output-dir 'data/${nextclade_dataset_name}'
    nextclade run -j $task.cpus --input-dataset data/${nextclade_dataset_name} --output-fasta ${name}.aligned.fasta --output-tsv tmp.tsv ${consensus}
    cat tmp.tsv | tr -d "\r" > ${name}_clade.tsv

    used_nextclade_version=\$nextclade_version_curr
    used_nextcladedataset_tag=\$(grep -Po '"tag":.*' data/${nextclade_dataset_name}/pathogen.json | cut -d' ' -f 2 | tr -d '"' | tr -d ',')
    used_nextcladedataset_info="${nextclade_dataset_name}, \$used_nextcladedataset_tag"
    """
    stub:
    """
    touch ${name}_clade.tsv ${name}.aligned.fasta
    used_nextclade_version=42
    used_nextcladedataset_info=42
    """
}
