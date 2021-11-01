include { filter_variants_hard; consensus_ambiguous; create_low_coverage_mask } from '../modules/bcftools'
include { adapt_consensus_header; consensus_masked } from '../modules/utils'

workflow generate_consensus{
    take:
        vcf
        reference
        bam

    main:
        filter_variants_hard(vcf)

        if (params.cns_gt_adjust > 0) {
            adjust_gt(filter_variants_hard.out)
        }
        vcf = params.cns_gt_adjust ? adjust_gt.out : filter_variants_hard.out

        adjust_del(vcf)
        consensus_ambiguous(adjust_del.out, reference, create_low_coverage_mask(bam).out) 
            adapt_consensus_header |
            consensus_masked

    emit:
        hard_filtered_variants = filter_variants.out
        consensus_ambiguous = adapt_consensus_header.out
        consensus_masked = consensus_masked.out
}