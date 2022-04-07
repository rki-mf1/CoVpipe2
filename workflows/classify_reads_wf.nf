include { kraken; filter_virus_reads } from '../modules/kraken'
include { krona; krona_taxonomy_update } from '../modules/krona'

workflow classify_reads {
    take:
        reads
        db
    main:
        kraken(reads, db)
        filter_virus_reads(kraken.out.fastq)

        // this has to be done once if the engine is conda/mamba; singularity works without it
        if ( workflow.profile.contains('mamba') ||  workflow.profile.contains('conda') ){
            krona_taxonomy_update()
            krona_tax_status = krona_taxonomy_update.out
        } else { krona_tax_status = 'not updated' }

        krona(kraken.out.kraken_report, krona_tax_status)
    emit:
        reads = filter_virus_reads.out.fastq
        report = kraken.out.kraken_report
}