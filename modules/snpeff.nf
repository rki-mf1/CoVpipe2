process snpeff {
    label 'snpeff'

    // publishDir "${params.output}/${params.variant_calling_dir}/${name}/", mode: params.publish_dir_mode

    input:
        tuple val(name), path(vcf)
        val(dataset_id)

    output:
        tuple val(name), path("${vcf.baseName}.annotation.html"), emit: html
        tuple val(name), path("${vcf.baseName}.annotation.covered.af.vcf"), emit: vcf

    script:
    """
    snpEff ann \
        -noLog \
        -stats ${vcf.baseName}.annotation.html \
        ${dataset_id} \
        ${vcf} 1> ${vcf.baseName}.annotation.covered.af.vcf
    """
    stub:
    """
    touch ${vcf.baseName}.annotation.html ${vcf.baseName}.annotation.covered.af.vcf
    """
}