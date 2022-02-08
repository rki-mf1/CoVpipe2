process freebayes {
    label 'freebayes' 

    // publishDir "${params.output}/${params.variant_calling_dir}/${name}/", mode: params.publish_dir_mode

    input:
    path(reference_fasta)
    path(reference_index)
    tuple val(name), path(bam), path(bai)

    output:
    tuple val(name), path("${bam.simpleName}.vcf")

    script:
    """
    freebayes -f ${reference_fasta} --min-alternate-count ${params.vcount} --min-alternate-fraction ${params.frac} --min-coverage ${params.cov} --pooled-continuous --haplotype-length -1 ${bam} | bcftools norm -f ${reference_fasta} -o ${bam.simpleName}.vcf
    """
    stub:
    """
    touch ${bam.simpleName}.vcf
    """
}