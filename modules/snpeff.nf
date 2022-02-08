process snpeff {
    label 'snpeff'

    // publishDir "${params.output}/${params.variant_calling_dir}/${name}/", mode: params.publish_dir_mode

    input:
        tuple val(name), path(vcf)
        path(reference)

    output:
        tuple val(name), path("${vcf.baseName}.annotation.html"), emit: html
        tuple val(name), path("${vcf.baseName}.annotation.covered.af.vcf"), emit: vcf

    script:
    """
    # get genome name   
    genome_name=\$(head -n1 ${reference} | cut -f1 -d' ' | sed 's/>//')

    snpEff ann \
        -noLog \
        -stats ${vcf.baseName}.annotation.html \
        \$genome_name \
        ${vcf} 1> ${vcf.baseName}.annotation.covered.af.vcf
    """
    stub:
    """
    touch ${vcf.baseName}.annotation.html ${vcf.baseName}.annotation.covered.af.vcf
    """
}