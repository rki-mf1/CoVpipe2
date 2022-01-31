include { kraken; filter_virus_reads; download_kraken_db } from '../modules/kraken'

workflow classify_reads {
    take:
        reads
    main:

        preload = file("${params.databases}/kraken2/kraken.tar.gz")
        if ( preload.exits() ) {
            db = preload
        } else {
            download_kraken_db()
            db = download_kraken_db.out
        }

        filter_virus_reads(kraken(reads, db).fastq)
    emit:
        reads = filter_virus_reads.out.fastq
        report = kraken.out.kraken_report
}