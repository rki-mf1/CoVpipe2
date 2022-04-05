include { index_bwa; bwa } from '../modules/bwa'
include { get_genomecov } from '../modules/bedtools'
include { index_bam; stats; get_fragment_size } from '../modules/samtools'

workflow mapping {
    take: 
        illumina_reads
        reference_fasta
    main:

        index_bwa(reference_fasta)
        bwa(illumina_reads, index_bwa.out) \
            | (index_bam & get_genomecov & stats & get_fragment_size)
        stats_output = stats.out.stats_small

    emit:
        bam_bai= index_bam.out
        coverage = get_genomecov.out.tsv
        mapping_stats = stats_output
        fragment_size = get_fragment_size.out
}
