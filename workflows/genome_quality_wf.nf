
include { president } from '../modules/president'

workflow genome_quality {
    take:
        fasta
        reference
    main:
        president(fasta, reference)
    emit:
        report = president.out.report
        valid = president.out.valid
        invalid = president.out.invalid
}