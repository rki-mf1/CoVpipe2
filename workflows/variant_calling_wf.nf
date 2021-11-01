include { freebayes } from '../modules/freebayes'

workflow variant_calling {
    take: 
        reference
        reference_fai
        bam_bai

    main:
        freebayes(reference, reference_fai, bam_bai)
        
    emit:
        vcf = freebayes.out
}
