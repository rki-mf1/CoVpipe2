include { freebayes }      from '../modules/freebayes'
include { bgzip_compress } from '../modules/utils'         addParams ( publish_dir: "${params.output}/${params.variant_calling_dir}/" )
include { index_vcf }      from '../modules/bcftools'      addParams ( publish_dir: "${params.output}/${params.variant_calling_dir}/" )

workflow variant_calling {
    take: 
        reference
        reference_fai
        bam_bai

    main:
        freebayes(reference, reference_fai, bam_bai) |
            bgzip_compress |
            index_vcf
        
    emit:
        vcf = freebayes.out
}
