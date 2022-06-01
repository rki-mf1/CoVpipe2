process nextclade {
    label 'nextclade'
    container = params.nextclade_docker
    if ( workflow.profile.contains('conda') ||  workflow.profile.contains('mamba') ) { conda = params.nextclade_conda }
    publishDir "${params.output}/${params.linage_dir}/${name}/", mode: params.publish_dir_mode

    input:
    tuple val(name), path(consensus)
    val(nextclade_dataset_name)


    output:
    tuple val(name), path("${name}_clade.tsv"), emit: results
    env(used_nextclade_version), emit: version
    env(used_nextcladedataset_version), emit: dataset_version

    script:
    """
    nextclade_version_curr=\$(nextclade --version)
    nextclade dataset get --name ${nextclade_dataset_name} --output-dir 'data/${nextclade_dataset_name}'
    nextclade run --input-fasta ${consensus} --input-dataset data/${nextclade_dataset_name} --output-tsv tmp.tsv
    cat tmp.tsv | tr -d "\r" > ${name}_clade.tsv

    used_nextclade_version=\$nextclade_version_curr
    used_nextcladedataset_version=\$(nextclade dataset list --name '${nextclade_dataset_name}' | grep 'Tag' | awk '{print \$3}')
    """
    stub:
    """
    touch ${name}_clade.tsv
    used_nextclade_version=42
    used_nextcladedataset_version=42
    """
}
