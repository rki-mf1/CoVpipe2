include { kraken; filter_virus_reads } from '../modules/kraken'
include { krona; krona_taxonomy_update } from '../modules/krona'
include { lcs_setup; lcs_sc2; lcs_plot } from '../modules/lcs'

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

        // calculate mixed/pooled samples using LCS, https://github.com/rvalieris/LCS
        if (params.read_linage) {
            lcs_setup()

            // lcs_sc2(filter_virus_reads.out.fastq)

            // lcs_results = lcs_sc2.out.map {it -> it[1]}.collectFile(name: 'lcs_results.tsv', skip: 1, keepHeader: true, storeDir: "${params.output}/${params.read_dir}/")

            // lcs_plot(lcs_results, params.lcs_cutoff)
        } else {
            lcs_output = Channel.empty()
        }

    emit:
        reads = filter_virus_reads.out.fastq
        report = kraken.out.kraken_report
}