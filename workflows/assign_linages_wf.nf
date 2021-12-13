include { pangolin ; update_pangolin } from '../modules/pangolin'

workflow assign_linages {
    take:
        fasta
    main:
        update_pangolin()
        pangolin(fasta, update_pangolin.out)
    emit:
        pangolin.out
}