process nextclade {
    label 'nextclade'
    container = params.nextclade_docker
    publishDir "${params.output}/${params.variant_calling_dir}/${name}/", mode: params.publish_dir_mode

    input:
    tuple val(name), path(consensus)
    val(nextclade_dataset_name)
    val(nextclade_version)

    output:
    tuple val(name), path("${name}_clade.tsv"), emit: results
    env(used_nextclade_version), emit: version
    env(used_nextcladedataset_version), emit: dataset_version

    script:
    """
    nextclade_version_curr=\$(nextclade --version)
    if [ "${nextclade_version}" != '' ]; then
        if [ "\$nextclade_version_curr" != "${nextclade_version}" ]; then
            echo "Something wrong in the nextclade update process."
            exit 1;
        fi
    fi
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

process update_nextclade_conda_env {
    // execute this locally - would most likely fail on custer systems
    label 'nextclade'
    executor 'local'
    cpus 1
    cache false

    output:
    env(nextclade_version)

    script:
    conda_mode = workflow.profile.contains('mamba') ? 'mamba' : 'conda'
    """
    ${conda_mode} update nextclade
    nextclade_version=\$(nextclade --version)
    """
    stub:
    """
    nextclade_version=42
    """
}