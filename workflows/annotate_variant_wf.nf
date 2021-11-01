include { snpeff } from '../modules/snpeff'

workflow annotate_variant {    
    take:
        vcf
        reference

    main:
        snpeff(vcf, reference)

    emit:
        html = snpeff.out.html
        vcf = snpeff.out.vcf
}