include { pangolin } from '../modules/pangolin'

workflow assign_linages {
    take:
        fasta
    main:
        pangolin(fasta)
    emit:
        pangolin.out
}