process bamclipper {
    label 'bamclipper'

    publishDir "${params.output}/${params.mapping_dir}/${name}/clipped", mode: params.publish_dir_mode

    input:
    tuple val(name), path(bam), path(bam_idx), path(primer_bed)

    output:
    tuple val(name), path("${bam.simpleName}.primerclipped.bam"), path("${bam.simpleName}.primerclipped.bam.bai")

    script:
    additional_parameters = params.bamclipper_additional_parameters ? params.bamclipper_additional_parameters : ''
    """
    bamclipper.sh -b ${bam} -p ${primer_bed} -n ${task.cpus} ${additional_parameters}
    """
    stub:
    """
    touch ${bam.simpleName}.primerclipped.bam ${bam.simpleName}.primerclipped.bam.bai
    """
}