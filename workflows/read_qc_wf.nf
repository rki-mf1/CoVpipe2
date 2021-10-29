include { fastp } from '../modules/fastp'

workflow read_qc {
    take: 
        illumina_reads
        adapter_fasta
    main:
        fastp(illumina_reads, adapter_fasta)
        reads_trimmed = fastp.out.reads
        fastp_json = fastp.out.json
    emit: 
        reads_trimmed
        fastp_json
}