include { index_bwa; bwa } from '../modules/bwa'
include { get_genomecov } from '../modules/bedtools'
include { index_bam; idxstats; flagstat; get_fragment_size } from '../modules/samtools'

workflow mapping {
    take: 
        illumina_reads
        reference_fasta
    main:

        index_bwa(reference_fasta)
        bwa(illumina_reads, index_bwa.out) \
            | (index_bam & get_genomecov & idxstats & flagstat & get_fragment_size)
        flagstat_output = flagstat.out
        idxstats_output = idxstats.out

    emit:
        bam_bai= index_bam.out
        coverage = get_genomecov.out.tsv
        flagstat = flagstat_output.flagstat
        flagstat_csv = flagstat_output.csv
        idxstats = idxstats_output
        fragment_size = get_fragment_size.out
}
