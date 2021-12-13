include { isec_vois } from '../modules/bcftools'
include { make_voi_table } from '../modules/utils'

workflow inspect_vois {
    take:
        vois
        vars
        low_coverage
    main:
        isec_vois(vois, vars.join(low_coverage)) \
            | make_voi_table
    emit:
        make_voi_table.out
}