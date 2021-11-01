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
}

process filter_variants_hard {
    label 'bcftools'

    input:
    tuple val(name), path(vcf)

    output:
    tuple val(name), path(vcf)

    script:
    """
    bcftools filter -e \
                    "INFO/MQM < ${params.var_mqm} | INFO/SAP > ${params.var_sap} | QUAL < ${params.var_qual}" \
                    -o ${name}.filtered.vcf -O z ${vcf}
    """
}

process create_low_coverage_mask {
    label 'bcftools'

    input:
    tuple val(name), file(bam)

    output:
    tuple val(name), file("${name}.lowcov.bed")

    script:
    """
    bedtools genomecov -bga -ibam ${bam} | awk '\$4 < ${params.cns_min_cov}' | bedtools merge  > ${name}.lowcov.bed
    """
}


process consensus_ambiguous {
    label 'bcftools'

    input:
    tuple val(name), path(vcf)
    path(reference)
    tuple val(name), path(mask_bed)


    output:
    path("${vcf.baseName}.iupac_consensus.tmp")

    script:
    """
    bcftools consensus \
                    -I \
                    -o ${vcf.baseName}.iupac_consensus.tmp \
                    -f ${reference} \
                    -m ${mask_bed} \
                    --sample ${vcf.baseName} \
                    ${vcf}
    """
}
