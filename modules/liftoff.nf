process liftoff {
    label 'liftoff'

    publishDir "${params.output}/${params.consensus_dir}/${name}", mode: params.publish_dir_mode

    input:
    tuple val(name), path(consensus)
    path(reference_genome)
    path(reference_annotation)

    output:
    path("${consensus.baseName}.gff"), emit: gff
    path("${consensus.baseName}_unmapped_features.txt"), emit: txt
    
    script:
    """
    liftoff -o ${consensus.baseName}.gff \
            -u ${consensus.baseName}_unmapped_features.txt -g ${reference_annotation} \
            ${consensus} ${reference_genome}
    """
    stub:
    """
    touch ${consensus.baseName}.gff ${consensus.baseName}_unmapped_features.txt
    """
} 