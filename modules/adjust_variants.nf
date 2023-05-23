process adjust_gt {
    label 'python'
    label 'smallTask'

    input:
    tuple val(name), path(vcf)
    val(cns_gt_adjust)

    output:
    tuple val(name), path("${name}.filtered.gt_adjust.vcf")

    script:
    """
    adjust_gt.py --vf ${cns_gt_adjust} --gz ${vcf} -o ${name}.filtered.gt_adjust.vcf
    """
    stub:
    """
    touch ${name}.filtered.gt_adjust.vcf
    """
}

process filter_indels {
    label 'python'
    label 'smallTask'

    input:
    tuple val(name), path(vcf)
    val(indel_threshold)

    output:
    tuple val(name), path("${name}.filtered.gt_adjust.filtered_indels.vcf")

    script:
    """
    filter_indels.py --vf ${indel_threshold} --gz ${vcf} -o ${name}.filtered.gt_adjust.filtered_indels.vcf
    """
    stub:
    """
    touch ${name}.filtered.gt_adjust.filtered_indels.vcf
    """
}
