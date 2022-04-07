include { pangolin ; update_pangolin_conda_env } from '../modules/pangolin'

workflow assign_linages {
    take:
        fasta
    main:
        if ( ( workflow.profile.contains('conda') || workflow.profile.contains('mamba') || workflow.profile.contains('standard') ) && params.update) {
            update_pangolin_conda_env()
            version = update_pangolin_conda_env.out
        } else {
            version = ''
        }
        pangolin(fasta, version)
    emit:
        report = pangolin.out.report
}