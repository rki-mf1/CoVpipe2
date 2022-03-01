include { pangolin ; update_pangolin } from '../modules/pangolin'

workflow assign_linages {
    take:
        fasta
    main:
        if (params.update_pangolin) {
            update_pangolin()
            version = update_pangolin.out
        } else {
            version = ''
        }
        pangolin(fasta, version)
    emit:
        report = pangolin.out.report
        version = pangolin.out.version
        scorpio_version = pangolin.out.scorpio_version
        scorpio_constellations_version = pangolin.out.scorpio_constellations_version
}