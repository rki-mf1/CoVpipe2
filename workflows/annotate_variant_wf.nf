include { snpeff }                                    from '../modules/snpeff'
include { nextclade }     from '../modules/nextclade'
include { bgzip_compress; replace_in_file }           from '../modules/utils'         addParams ( publish_dir: "${params.output}/${params.variant_calling_dir}/" )
include { index_vcf }                                 from '../modules/bcftools'      addParams ( publish_dir: "${params.output}/${params.variant_calling_dir}/" )

workflow annotate_variant {    
    take:
        vcf
        consensus
        nextclade_dataset_name
        snpeff_dataset_id

    main:
        replace_in_file(vcf, 'MN908947.3', snpeff_dataset_id) // really sars-cov-2 specific
        snpeff(replace_in_file.out, snpeff_dataset_id)

        bgzip_compress(snpeff.out.vcf) \
            | index_vcf

        nextclade(consensus, nextclade_dataset_name)

    emit:
        html = snpeff.out.html
        vcf = snpeff.out.vcf
        nextclade_results = nextclade.out.results
        nextclade_version = nextclade.out.version
        nextclade_dataset_version = nextclade.out.dataset_version
}
