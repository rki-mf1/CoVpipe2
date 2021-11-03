include { liftoff } from '../modules/liftoff'

workflow annotate_consensus {
    take:
        consensus
        reference_genome
        reference_annotation
    main:
        liftoff(consensus, reference_genome, reference_annotation)
    emit:
        liftoff.out.gff
}