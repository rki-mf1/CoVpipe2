process desh_qc {
    label 'bcftools'

    input:
    tuple val(name), path(consensus)
    tuple val(name), path(coverage)
    val(min_cov)

    output: 
    tuple val(name), path("${name}_desh.tsv")

    script:
    """
    echo 'sample,#n,#iupac,#lowcov' > ${name}_desh.tsv
    N=\$(tail -n +2 ${consensus} | tr -cd "N" | wc -c)
    IUPAC=\$(tail -n +2 ${consensus} | tr -cd "RYSWKMBDHVN" | wc -c)
    COVTH=\$(awk -F "\t" '{{COV=\$3; if (COV < ${min_cov}) {{ print COV }} }}' ${coverage} | wc -l)
    echo ${name},\$N,\$IUPAC,\$COVTH >> ${name}_desh.tsv
    """
}

// process fastp_table {
//     label 'r'
//     label 'smallTask'

//     input:
    

//     output:

//     script:
//     """
//     """
// }

// process kraken_table {
    // label 'r'
    // label 'smallTask'

    // input:

    // output:

    // script:
    // """
    // """
// }

process flagstat_table {
    label 'r'
    label 'smallTask'

    input:
    path(flagstat_csv)

    output:
    path("mapping_stats.csv")

    script:
    name_list = flagstat_csv.collect{ "\"${it.getBaseName()}\"" }.join(",")
    file_list = flagstat_csv.collect{ "\"${it}\"" }.join(",")
    """
    #!/usr/bin/env Rscript

    library("data.table")
    library("plyr")

    f.list <- c(${file_list})
    names(f.list) <- c(${name_list})

    df.bamstat.data <- ldply(f.list, fread, sep = ';')
    colnames(df.bamstat.data) <- c("sample", "count", "unknown", "description")

    df.output <- data.frame("sample" = unique(df.bamstat.data\$sample),
                            "input" = df.bamstat.data\$count[grepl("in total", df.bamstat.data\$description)],
                            "mapped" = df.bamstat.data\$count[grepl("properly paired", df.bamstat.data\$description)])

    df.output\$mapping.rate <- df.output\$mapped / df.output\$input

    # save table as csv for later use
    write.csv(  x=df.output, 
                row.names = FALSE, 
                file = file.path("mapping_stats.csv")
    )

    #df.output\$input <- f.color_bar("lightgreen")(df.output\$input)
    #df.output\$mapped <- f.color_bar("lightgreen")(df.output\$mapped)
    """
}

// process fragment_size_table {
    // label 'r'
    // label 'smallTask'

    // input:

    // output:

    // script:
    // """
    // """
// }

// process coverage_table {
    // label 'r'
    // label 'smallTask'

    // input:

    // output:

    // script:
    // """
    // """
// }

process multiqc_report {
    label 'multiqc'
    label 'smallTask'

    input:
    path(fastp)
    path(kraken)
    path(flagstat)
    path(pangolin)

    output:
    path ("multiqc_report.html")

    script:
    """
    multiqc .
    """
}

// process rmarkdown_report {
    // label 'r'
    // label 'smallTask'

    // input:

    // output:

    // script:
    // """
    // """
// }