include { bgzip_compress as bgzip_compress1; bgzip_compress as bgzip_compress2 } from '../modules/utils' addParams ( publish_dir: "${params.output}/${params.variant_calling_dir}/" )
include { adapt_consensus_header; mask_iupac as consensus_masked } from '../modules/utils' addParams ( publish_dir: "${params.output}/${params.consensus_dir}/" )
include { filter_variants_hard; index_vcf; vcf_filter_dels } from '../modules/bcftools' addParams ( publish_dir: "${params.output}/${params.variant_calling_dir}/" )
include { consensus_ambiguous } from '../modules/bcftools' addParams ( publish_dir: "${params.output}/${params.consensus_dir}/" )
include { adjust_gt; filter_indels } from '../modules/adjust_variants'
include { create_low_coverage_no_del_mask } from '../modules/bedtools' addParams ( publish_dir: "${params.output}/${params.consensus_dir}/" )

workflow generate_consensus{
    take:
        vcf
        reference
        bam

    main:
        filter_variants_hard(vcf)

        if (params.cns_gt_adjust > 0) {
            adjust_gt(filter_variants_hard.out, params.cns_gt_adjust) \
                | bgzip_compress1
                | set { gt_adjusted_vcf }
        }
        vcf = params.cns_gt_adjust > 0 ? gt_adjusted_vcf : filter_variants_hard.out
        
        if (params.cns_indel_filter > 0) {
            filter_indels(vcf, params.cns_indel_filter) \
                | bgzip_compress2
                | set { indel_filtered_vcf }
        }

        vcf = params.cns_indel_filter > 0 ? indel_filtered_vcf : vcf

        vcf_filter_dels (vcf) \
            | join(bam)
            | create_low_coverage_no_del_mask
            | set { low_cov_mask }
            
        index_vcf(vcf) \
            | join(low_cov_mask)
            | set { vfc_and_low_cov_mask }
        
        consensus_ambiguous(vfc_and_low_cov_mask, reference) \
            | adapt_consensus_header \
            | consensus_masked

    emit:
        hard_filtered_variants = filter_variants_hard.out
        consensus_ambiguous = adapt_consensus_header.out
        consensus_masked = consensus_masked.out
        low_coverage_bed = create_low_coverage_no_del_mask.out
}