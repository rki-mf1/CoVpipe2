
include { president } from '../modules/president'

workflow genome_quality {
    take:
        fasta
        reference
    main:
        president(fasta, reference)
        president_result = president.out.report
            | map {it -> it[1]} \
            | collectFile(name: 'president_results.tsv', skip: 1, keepHeader: true)
    emit:
        president_result
}