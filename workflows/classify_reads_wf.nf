include { kraken; filter_virus_reads } from '../modules/kraken'

workflow classify_reads {
    take:
        reads
        db
    main:
        filter_virus_reads(kraken(reads, db).fastq)
    emit:
        reads = filter_virus_reads.out.fastq
        report = kraken.out.kraken_report
}