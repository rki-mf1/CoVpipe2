process nextclade {
    label 'nextclade'

    // publishDir "${params.output}/${params.lineagedir}/${name}/", mode: 'copy', pattern: "${name}_clade.tsv"
    input:
    tuple val(name), path(consensus)
    val(nextclade_version)

    output:
    tuple val(name), path("${name}_clade.tsv")

    script:
    """
    nextclade_version_curr=\$(nextclade --version)
    if [ "${nextclade_version}" != '' ]; then
        if [ "\$nextclade_version_curr" != "${nextclade_version}" ]; then
            echo "Something wrong in the nextclade update process."
            exit 1;
        fi
    fi
    nextclade dataset get --name 'sars-cov-2' --output-dir 'data/sars-cov-2'
    nextclade run --input-fasta ${consensus} --input-dataset data/sars-cov-2 --output-tsv tmp.tsv
    cat tmp.tsv | tr -d "\r" > ${name}_clade.tsv
    """
    stub:
    """
    touch ${name}_clade.tsv
    """
}

process update_nextclade {
    // execute this locally - would most likely fail on custer systems
    label 'nextclade'
    executor 'local'
    cpus 1
    cache false

    output:
    env(nextclade_version)

    script:
    """
    conda update nextclade
    nextclade_version=\$(nextclade --version)
    """
}