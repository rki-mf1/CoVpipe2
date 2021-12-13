
include { president } from '../modules/president'

workflow genome_quality {
    take:
        fasta
        reference
    main:
        president(fasta, reference)
    emit:
        president.out.report
}