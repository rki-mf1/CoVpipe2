include { index_bwa; bwa } from '../modules/bwa'
include { index_hisat2; hisat2 } from '../modules/hisat2'
include { get_genomecov } from '../modules/bedtools'
include { index_bam; flagstat; get_fragment_size } from '../modules/samtools'

workflow mapping {
    take: 
        illumina_reads
        reference_fasta
    main:

        index_bwa(reference_fasta)
        bwa(illumina_reads, index_bwa.out) 

        index_hisat2(reference_fasta)
        hisat2(illumina_reads, index_hisat2.out)

        hisat2.out.bam.concat(bwa.out.bam) \
            | (index_bam & get_genomecov & flagstat & get_fragment_size)
        flagstat_output = flagstat.out
    emit:
        bam_bai= index_bam.out
        coverage = get_genomecov.out.tsv
        flagstat = flagstat_output.flagstat
        flagstat_csv = flagstat_output.csv
        fragment_size = get_fragment_size.out
}
