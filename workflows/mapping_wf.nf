include { index_bwa; bwa } from '../modules/bwa'
include { get_genomecov } from '../modules/bedtools'
include { index_bam; stats; get_fragment_size; filter_isize_bam } from '../modules/samtools'

workflow mapping {
    take: 
        illumina_reads
        reference_fasta
        insert_size_filter
    main:

        index_bwa(reference_fasta)
        bwa(illumina_reads, index_bwa.out)
        
        if (insert_size_filter != false) {
            filter_isize_bam(bwa.out, insert_size_filter)
        }
        mapping = insert_size_filter == false ? bwa.out.bam : filter_isize_bam.out.bam

        index_bam(mapping)
        get_genomecov(mapping)
        stats(mapping)
        get_fragment_size(mapping)

    emit:
        bam_bai = index_bam.out
        coverage = get_genomecov.out.tsv
        mapping_stats = stats.out.stats_small
        fragment_size = get_fragment_size.out
}
