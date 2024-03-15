
include { president } from '../modules/president'

workflow genome_quality {
    take:
        fasta
        reference
        seq_threshold
        n_threshold
    main:
        president(fasta, reference, seq_threshold, n_threshold)
    emit:
        valid = president.out.valid
        invalid = president.out.invalid
}