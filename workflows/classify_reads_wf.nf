include { kraken; filter_virus_reads } from '../modules/kraken'
include { krona } from '../modules/krona'


workflow classify_reads {
    take:
        reads
        db
    main:
        kraken(reads, db)
        filter_virus_reads(kraken.out.fastq)
        krona(kraken.out.kraken_report)
    emit:
        reads = filter_virus_reads.out.fastq
        un_reads = kraken.out.un_fastq
        report = kraken.out.kraken_report
}