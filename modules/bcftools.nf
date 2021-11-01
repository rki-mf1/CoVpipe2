process index_vcf {
    label 'bcftools'

    publishDir "${params.publish_dir}/${name}", mode: params.publish_dir_mode

    input:
    tuple val(name), path(vcf)

    output:
    tuple val(name), path(vcf), path("${vcf}.csi")

    script:
    """
    bcftools index -f ${vcf}
    """
}