process multiqc {
    conda 'bioconda:multiqc=1.12'
    publishDir "${params.output}/${params.report_dir}/", mode: params.publish_dir_mode

    input:
        path(fastqc)
        path(krona)
        path(mapping)
    output:
        path("multiqc_report.html")

    script:
    """
    multiqc .
    """
}