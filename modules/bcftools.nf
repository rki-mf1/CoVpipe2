process index_vcf {
    label 'bcftools'

    publishDir "${params.publish_dir}/${name}", mode: params.publish_dir_mode

    input:
    tuple val(name), path(vcf)

    output:
    tuple val(name), path(vcf), path("${vcf}.csi")

    script:
    """
    bcftools index -f ${vcf}
    """
    stub:
    """
    touch "${vcf}.csi"
    """
}

process filter_variants_hard {
    label 'bcftools'

    publishDir "${params.output}/${params.consensus_dir}/${name}", mode: params.publish_dir_mode

    input:
    tuple val(name), path(vcf)

    output:
    tuple val(name), path("${name}.filtered.vcf.gz")

    script:
    var_sap_filter = params.var_sap == '-1' ? '' : "| INFO/SAP > ${params.var_sap}"
    """
    bcftools filter -e \
                    "INFO/MQM < ${params.var_mqm} ${var_sap_filter} | QUAL < ${params.var_qual}" \
                    -o ${name}.filtered.vcf.gz -O z ${vcf}
    """
    stub:
    """
    touch ${name}.filtered.vcf.gz
    """
}

process vcf_filter_dels {
    label 'bcftools'

    input:
    tuple val(name), path(vcf)

    output:
    tuple val(name), path("${vcf}.all_dels")

    script:
    """
    bcftools filter -i "INFO/TYPE=='del'" ${vcf} > ${vcf}.all_dels
    """
    stub:
    """
    touch ${vcf}.all_dels
    """
}

process consensus_ambiguous {
    label 'bcftools'

    input:
    tuple val(name), path(vcf), path(csi), path(mask_bed)
    path(reference)

    output:
    tuple val(name), path("${vcf.baseName}.iupac_consensus.tmp")

    script:
    """
    bcftools consensus \
                    -I \
                    -o ${vcf.baseName}.iupac_consensus.tmp \
                    -f ${reference} \
                    -m ${mask_bed} \
                    --sample ${name} \
                    ${vcf}
    """
    stub:
    """
    touch ${vcf.baseName}.iupac_consensus.tmp
    """
}

process isec_vois {
    label 'bcftools'

    input:
    path(vois)
    tuple val(name), path(vars), path(vars_csi), path(low_cov)
    
    output:
    path(vois)
    tuple val(name), path("${name}.voi_exact.vcf")
    tuple val(name), path("${name}.voi_not_found.vcf")
    tuple val(name), path("${name}.voi_low_coverage.vcf")
    tuple val(name), path("${name}.voi_diff_voi.vcf")
    tuple val(name), path("${name}.voi_diff_sample.vcf")

    script:
    """
    bgzip -c ${vois} > ${vois}.gz
    bcftools index -f ${vars}
    bcftools index -f ${vois}.gz


    bcftools isec -c none -n=2 -w2 -o ${name}.voi_exact.vcf ${vars} ${vois}.gz # identical REF and ALT
    bcftools isec -c all -n~01 -w2 -o ${name}.voi_not_found.vcf ${vars} ${vois}.gz # voi not found, regardless ALT
    bcftools isec -T ${low_cov} -o ${name}.voi_low_coverage.vcf ${vois}.gz # low coverage vois

    bgzip -c ${name}.voi_exact.vcf > ${name}.voi_exact.vcf.gz
    bcftools index -f ${name}.voi_exact.vcf.gz
    bgzip -c ${name}.voi_not_found.vcf > ${name}.voi_not_found.vcf.gz
    bcftools index -f ${name}.voi_not_found.vcf.gz
    
    bcftools isec -n~100 -w1 -o ${name}.voi_diff_voi.vcf ${vois}.gz ${name}.voi_exact.vcf.gz ${name}.voi_not_found.vcf.gz # vio with different ALT in voi file

    bgzip -c ${name}.voi_diff_voi.vcf > ${name}.voi_diff_voi.vcf.gz
    bcftools index -f ${name}.voi_diff_voi.vcf.gz
    
    bcftools isec -c all -n=2 -w1 -o ${name}.voi_diff_sample.vcf ${vars} ${name}.voi_diff_voi.vcf.gz # voi with different ALT in vcf sample file
    """
    stub:
    """
    touch ${name}.voi_exact.vcf ${name}.voi_not_found.vcf ${name}.voi_low_coverage.vcf ${name}.voi_diff_voi.vcf ${name}.voi_diff_sample.vcf
    """
}