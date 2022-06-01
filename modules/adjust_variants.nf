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
