include { snpeff }         from '../modules/snpeff'
include { nextclade; update_nextclade }      from '../modules/nextclade'
include { bgzip_compress } from '../modules/utils'         addParams ( publish_dir: "${params.output}/${params.variant_calling_dir}/" )
include { index_vcf }      from '../modules/bcftools'      addParams ( publish_dir: "${params.output}/${params.variant_calling_dir}/" )

workflow annotate_variant {    
    take:
        vcf
        consensus
        reference

    main:
        snpeff(vcf, reference)

        bgzip_compress(snpeff.out.vcf) \
            | index_vcf

        if (params.update_nextclade) {
            update_nextclade()
            version = update_nextclade.out
        } else {
            version = ''
        }
        nextclade(consensus, version)

    emit:
        html = snpeff.out.html
        vcf = snpeff.out.vcf
        results = nextclade.out.results
        nextclade_version = nextclade.out.version
        nextclade_dataset_version = nextclade.out.dataset_version
}
