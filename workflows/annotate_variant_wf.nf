include { snpeff }         from '../modules/snpeff'
include { bgzip_compress } from '../modules/utils'         addParams ( publish_dir: "${params.output}/${params.variant_calling_dir}/" )
include { index_vcf }      from '../modules/bcftools'      addParams ( publish_dir: "${params.output}/${params.variant_calling_dir}/" )

workflow annotate_variant {    
    take:
        vcf
        reference

    main:
        snpeff(vcf, reference)

        bgzip_compress(snpeff.out.vcf) |
            index_vcf

    emit:
        html = snpeff.out.html
        vcf = snpeff.out.vcf
}