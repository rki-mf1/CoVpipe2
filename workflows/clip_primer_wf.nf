include { bamclipper } from '../modules/bamclipper'

workflow clip_primer {
    take:
        bam_bai
        primer_bed

    main:
        bamclipper(bam_bai, primer_bed)

    emit:
        bam_bai = bamclipper.out
}