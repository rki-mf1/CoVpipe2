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

process adjust_del {
    label 'python'
    label 'smallTask'
    
    input:
    tuple val(name), path(vcf)

    output:
    tuple val(name), path("${vcf.baseName}.del_adjusted.vcf")

    script:
    """
    adjust_gt.py --gz ${vcf} -o ${vcf.baseName}.del_adjusted.vcf
    """
    stub:
    """
    touch ${vcf.baseName}.del_adjusted.vcf
    """
}