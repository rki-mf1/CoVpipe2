process adjust_gt {
    label 'python'
    label 'smallTask'

    input:
    tuple val(name), path(vcf)

    output:
    tuple val(name), path("${name}.filtered.gt_adjust.vcf")

    script:
    """
    python ../bin/adjust_gt.py --vf ${params.cns_gt_adjust} --gz ${vcf} -o ${name}.filtered.gt_adjust.vcf
    """
}

process adjust_del {
    label 'python'
    label 'smallTask'
    
    input:
    tuple val(name), path(vcf)

    output:
    tuple val(name), path()

    script:
    """
    python ../bin/adjust_gt.py --gz ${vcf} -o ${name}.filtered.del_adjusted.vcf
    """
}