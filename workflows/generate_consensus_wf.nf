include { bgzip_compress as bgzip_compress_1; bgzip_compress as bgzip_compress_2; bgzip_compress as bgzip_compress_3; adapt_consensus_header; mask_iupac as consensus_masked } from '../modules/utils' addParams ( publish_dir: "${params.output}/${params.consensus_dir}/" )
include { filter_variants_hard; consensus_ambiguous; index_vcf } from '../modules/bcftools' addParams ( publish_dir: "${params.output}/${params.consensus_dir}/" )
include { adjust_gt; adjust_del } from '../modules/adjust_variants'
include { create_low_coverage_mask } from '../modules/bedtools' addParams ( publish_dir: "${params.output}/${params.consensus_dir}/" )

workflow generate_consensus{
    take:
        vcf
        reference
        bam

    main:
        filter_variants_hard(vcf)

        if (params.cns_gt_adjust > 0) {
            adjust_gt(filter_variants_hard.out, params.cns_gt_adjust) \
                | bgzip_compress_2
                | set { gt_adjusted_vcf }
        }
        vcf = params.cns_gt_adjust > 0 ? gt_adjusted_vcf : filter_variants_hard.out

        adjust_del(vcf) \
            | bgzip_compress_3 \
            | index_vcf \
            | set { del_adjusted_vcf }
        
        consensus_ambiguous(del_adjusted_vcf.join(create_low_coverage_mask(bam)), reference) \
            | adapt_consensus_header \
            | consensus_masked

    emit:
        hard_filtered_variants = filter_variants_hard.out
        consensus_ambiguous = adapt_consensus_header.out
        consensus_masked = consensus_masked.out
        low_coverage_bed = create_low_coverage_mask.out
}