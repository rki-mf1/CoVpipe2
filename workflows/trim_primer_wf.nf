include { ptrimmer } from '../modules/ptrimmer'
include { compress_reads } from '../modules/utils'

workflow trim_primer {
    take:
        illumina_reads
        primer_set
    main:
        ptrimmer(illumina_reads, primer_set)
        compress_reads(ptrimmer.out.reads)
    emit:
        compress_reads.out
}