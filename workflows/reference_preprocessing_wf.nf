include { dos2unix } from '../modules/dos2unix'
include { index_fasta } from '../modules/samtools'

workflow reference_preprocessing {
    take: reference_fasta
    main:
        dos2unix(reference_fasta) \
            | index_fasta
    emit: 
        ref = dos2unix.out
        fai = index_fasta.out
}