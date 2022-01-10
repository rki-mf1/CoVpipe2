include { pangolin ; update_pangolin } from '../modules/pangolin'

workflow assign_linages {
    take:
        fasta
    main:
        if (params.update_pangolin) {
            update_pangolin()
            pangolin_version = update_pangolin.out
        } else {
            pangolin_version = ''
        }
        pangolin(fasta, pangolin_version)
    emit:
        pangolin.out
}