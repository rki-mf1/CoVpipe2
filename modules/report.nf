process fastp_table {
    label 'r'

    input:
    path(json)

    output:
    path("read_stats.csv"), emit: stats
    path("read_stats_filter.csv"), emit: stats_filter

    script:
    name_list = json.collect{ "\"${it.getSimpleName()}\"" }.join(",")
    file_list = json.collect{ "\"${it}\"" }.join(",")
    """
    #!/usr/bin/env Rscript

    library("data.table")
    library("plyr")
    library("rjson")

    f.list <- c(${file_list})
    names(f.list) <- c(${name_list})

    l.trimming.data.json <- lapply(f.list, function(x){
        fromJSON(file = x)
    })

    df.trimming.data <- ldply(l.trimming.data.json, function(e){
        df.before <- as.data.frame(do.call(rbind, e\$summary\$before_filtering))
        colnames(df.before) <- c("before.trimming")
        df.after  <- as.data.frame(do.call(rbind, e\$summary\$after_filtering))
        colnames(df.after) <- c("after.trimming")

        df.output <- data.frame(feature = rownames(df.before),
                                before = df.before\$before.trimming,
                                after = df.after\$after.trimming)
        return(df.output)
    })

    # rename 1st column
    tmp <- colnames(df.trimming.data)
    tmp[1] <- c("sample")
    colnames(df.trimming.data) <- tmp

    df.filter.data <- ldply(l.trimming.data.json, function(e){
        
        df.output <- data.frame(passed_filter = e\$filtering_result\$passed_filter_reads,
                                low_qual = e\$filtering_result\$low_quality_reads,
                                high_N = e\$filtering_result\$too_many_N_reads,
                                low_complex = e\$filtering_result\$low_complexity_reads,
                                short = e\$filtering_result\$too_short_reads
                                )
        return(df.output)
    })

    # rename 1st column
    tmp <- colnames(df.filter.data)
    tmp[1] <- c("sample")
    colnames(df.filter.data) <- tmp

    df.summary <- data.frame(sample = unique(df.trimming.data\$sample),
                        reads.before.clip = df.trimming.data\$before[grepl("total_reads", df.trimming.data\$feature)],
                        reads.after.clip  = df.trimming.data\$after[grepl("total_reads", df.trimming.data\$feature)],
                        ratio.clip = df.trimming.data\$after[grepl("total_reads", df.trimming.data\$feature)] / 
                            df.trimming.data\$before[grepl("total_reads", df.trimming.data\$feature)],
                        q30.before.clip = df.trimming.data\$before[grepl("q30_rate", df.trimming.data\$feature)],
                        q30.after.clip  = df.trimming.data\$after[grepl("q30_rate", df.trimming.data\$feature)]
    )

    df.table <- df.summary

    # add coloured bar charts to table
    # df.table\$reads.before.clip <- f.color_bar("lightgreen")(df.table\$reads.before.clip)
    # df.table\$reads.after.clip <- f.color_bar("lightgreen")(df.table\$reads.after.clip)

    write.csv(  x=df.summary, 
            row.names = FALSE, 
            file = file.path("read_stats.csv")
        )
    write.csv(  x=df.filter.data, 
            row.names = FALSE, 
            file = file.path("read_stats_filter.csv")
        )
    """
    stub:
    """
    touch read_stats.csv read_stats_filter.csv
    """
}

process kraken_table {
    label 'r'

    input:
    path(report)
    val(tax_id)

    output:
    path("species_filtering.csv")

    script:
    name_list = report.collect{ "\"${it.getSimpleName()}\"" }.join(",")
    file_list = report.collect{ "\"${it}\"" }.join(",")

    """
    #!/usr/bin/env Rscript

    library("data.table")
    library("plyr")

    f.list <- c(${file_list})
    names(f.list) <- c(${name_list})

    df.kraken_output <- data.frame()

    dt.kraken_data <- ldply(f.list, fread)
    colnames(dt.kraken_data) <- c("sample", "read_ratio", "read_count", "read_count_specific", "rank", "ncbi_taxid", "sciname")
    # select unclassified and user supplied tax id
    df.kraken_output <- dt.kraken_data[dt.kraken_data\$ncbi_taxid %in% c(0, ${tax_id}, "9606"), c("sample", "read_ratio", "read_count", "ncbi_taxid", "sciname")]
    # make tables wide for absolute and relative read counts
    dt.ratio <- data.table::dcast(as.data.table(df.kraken_output), sample ~ sciname, value.var = c("read_ratio"))
    dt.count <- data.table::dcast(as.data.table(df.kraken_output), sample ~ sciname, value.var = c("read_count"))
    # join both tables for output
    df.kraken_output <- join(dt.ratio, dt.count, by = "sample")
    
    if (length(colnames(df.kraken_output)) == 3) {
        colnames(df.kraken_output) <- paste0(colnames(df.kraken_output), c("", " (ratio)", " (count)"))
    } else if (length(colnames(df.kraken_output)) == 5) {
        colnames(df.kraken_output) <- paste0(colnames(df.kraken_output), c("", " (ratio)", " (ratio)", " (count)", " (count)"))
    } else if (length(colnames(df.kraken_output)) == 7) {
        colnames(df.kraken_output) <- paste0(colnames(df.kraken_output), c("", " (ratio)", " (ratio)", " (ratio)", " (count)", " (count)", " (count)"))
    }

    # add rows for reproducible colour scaling
    if (length(colnames(df.kraken_output)) == 3) {
        df.tmp <- as.data.frame(matrix(rep(c(0,100),3), ncol = 3))
    } else if (length(colnames(df.kraken_output)) == 5) {
        df.tmp <- as.data.frame(matrix(rep(c(0,100),5), ncol = 5))
    } else if (length(colnames(df.kraken_output)) == 7) {
        df.tmp <- as.data.frame(matrix(rep(c(0,100),7), ncol = 7))
    }
    colnames(df.tmp) <- colnames(df.kraken_output)
    df.kraken_output <- rbind(df.kraken_output, df.tmp)
    
    # remove added lines for colour scaling
    df.kraken_output <- head(df.kraken_output, n = -2)

    # save table as csv for later use
    write.csv(  x=df.kraken_output, 
                row.names = FALSE, 
                file = file.path("species_filtering.csv")
    )
    """
    stub:
    """
    touch species_filtering.csv
    """
}

process flagstat_table {
    label 'r'

    input:
    path(flagstat_csv)

    output:
    path("mapping_stats.csv")

    script:
    name_list = flagstat_csv.collect{ "\"${it.getSimpleName()}\"" }.join(",")
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

    write.csv(  x=df.output, 
                row.names = FALSE, 
                file = file.path("mapping_stats.csv")
    )

    #df.output\$input <- f.color_bar("lightgreen")(df.output\$input)
    #df.output\$mapped <- f.color_bar("lightgreen")(df.output\$mapped)
    """
    stub:
    """
    touch mapping_stats.csv
    """
}

process fragment_size_table {
    label 'r'

    input:
    path(tsv)

    output:
    path("fragment_sizes.csv"), emit: size
    path("fragment_sizes_median.csv"), emit: median

    script:
    name_list = tsv.collect{ "\"${it.getSimpleName()}\"" }.join(",")
    file_list = tsv.collect{ "\"${it}\"" }.join(",")
    """
    #!/usr/bin/env Rscript

    library("data.table")
    library("plyr")

    f.list <- c(${file_list})
    names(f.list) <- c(${name_list})

    # read.csv produces a new column for each read file
    dt.fragsizes <- as.data.table(ldply(f.list, fread))
    colnames(dt.fragsizes) <- c("sample", "fragsize")

    dt.fragsizes\$fragsize.abs <- abs(dt.fragsizes\$fragsize)
    dt.fragsizes.median <- dt.fragsizes[, median(fragsize.abs), by = c("sample")]
    setnames(dt.fragsizes.median, "V1", "median.fragsize")
    dt.fragsizes.median\$median.fragsize <- as.numeric(as.character(dt.fragsizes.median\$median.fragsize))

    # add standard deviation
    dt.fragsizes.median\$sd.fragsize <- dt.fragsizes[, sd(fragsize.abs), by = c("sample")]\$V1

    # conditional formatting for fragsizes not between 90 and 110
    #dt.fragsizes.median\$median.fragsize <- ifelse(dt.fragsizes.median\$median.fragsize >= 110 | dt.fragsizes.median\$median.fragsize <= 90,
    #                                          cell_spec(dt.fragsizes.median\$median.fragsize, background = "orange", align = "right"),
    #                                          dt.fragsizes.median\$median.fragsize)

    write.csv(  x=dt.fragsizes.median, 
            row.names = FALSE, 
            file = file.path("fragment_sizes_median.csv")
    )
    write.csv(  x=dt.fragsizes, 
            row.names = FALSE, 
            file = file.path("fragment_sizes.csv")
    )
    """
    stub:
    """
    touch fragment_sizes.csv fragment_sizes_median.csv
    """
}

process coverage_table {
    label 'r'

    input:
    path(tsv)
    val(min_cov)

    output:
    path("coverage_table.csv"), emit: coverage_table
    path("positive_samples.csv"), emit: positive
    path("negative_samples.csv"), emit: negative
    path("coverage_samples.csv"), emit: sample_cov

    script:
    name_list = tsv.collect{ "\"${it.getSimpleName()}\"" }.join(",")
    file_list = tsv.collect{ "\"${it}\"" }.join(",")
    """
    #!/usr/bin/env Rscript

    library("data.table")
    library("plyr")
    
    f.list <- c(${file_list})
    names(f.list) <- c(${name_list})

    dt.coverage <- as.data.table(ldply(f.list, fread))
    colnames(dt.coverage) <- c("sample","chromosome", "position", "depth")

    # reduce amount of data points to be plotted
    dt.coverage[, bin:=rep(seq(1, ceiling(length(position) / 100)), each = 100, length.out = length(position)), by = "sample"]
    dt.coverage[, mid.bin:=seq(1,length(position)) %% 100, by = "sample"] # by samples is here also necessary
    dt.coverage[, mean.cov:=mean(depth), by=c("sample", "bin")]

    dt.output <- dt.coverage[, sum(depth > ${min_cov}), by = sample] #
    setnames(dt.output, "V1", "covered.bases")

    dt.output\$genome.length <- dt.coverage[,length(depth), by = sample]\$V1
    dt.output\$genome.coverage <- dt.output\$covered.bases / dt.output\$genome.length

    dt.output\$DP.median <- dt.coverage[, median(depth), by = sample]\$V1
    dt.output\$DP.mean <- dt.coverage[, mean(depth), by = sample]\$V1

    write.csv( x=dt.coverage, 
            row.names = FALSE, 
            file = file.path("coverage_table.csv"))
    # positve and negative tables as by consideration of cleanplex quality measures (95% reference genome coverage)
    write.csv( x=dt.output[dt.output\$genome.coverage >= 0.95], 
            row.names = FALSE, 
            file = file.path("positive_samples.csv"))
    write.csv( x=dt.output[dt.output\$genome.coverage < 0.95], 
            row.names = FALSE, 
            file = file.path("negative_samples.csv"))
    # complete table output
    write.csv( x=dt.output, 
            row.names = FALSE, 
            file = file.path("coverage_samples.csv"))
    """
    stub:
    """
    touch coverage_table.csv positive_samples.csv negative_samples.csv coverage_samples.csv
    """
}

process rmarkdown_report {
    // rmarkdown::render does not respect symlinks https://github.com/rstudio/rmarkdown/issues/1508
    label 'r'

    publishDir "${params.output}/${params.report_dir}/", mode: params.publish_dir_mode
    // stageInMode: 'copy' # would copy every input file, not so nice

    input:
    path(rmd)
    path(fastp_table_stats)
    path(fastp_table_stats_filter)
    path(kraken_table)
    path(flagstat_table)
    path(fragment_size_table)
    path(fragment_size_median_table)
    path(coverage_table)
    path(positive)
    path(negative)
    path(sample_cov)
    path(president_results)
    path(pangolin_results)
    val(pangolin_version)
    path(nextclade_results)
    path(vois_results)

    output:
    path("report.html")

    script:
    kraken_table_optional = kraken_table ? kraken_table : 'none'
    vois_results_optional = vois_results ? vois_results : 'none'
    run_id = params.run_id != '' ? params.run_id : 'none'
    pipeline_version = workflow.repository != null ? "$workflow.repository - $workflow.revision [$workflow.commitId]" : 'none'
    """
    cp -L ${rmd} report.Rmd
    Rscript -e "rmarkdown::render('report.Rmd', params=list(fastp_table_stats='${fastp_table_stats}', fastp_table_stats_filter='${fastp_table_stats_filter}', kraken_table='${kraken_table_optional}', flagstat_table='${flagstat_table}', fragment_size_table='${fragment_size_table}', fragment_size_median_table='${fragment_size_median_table}', coverage_table='${coverage_table}', positive='${positive}', negative='${negative}', sample_cov='${sample_cov}', president_results='${president_results}', pangolin_results='${pangolin_results}', pangolin_version='${pangolin_version}', nextclade_results='${nextclade_results}', vois_results='${vois_results_optional}', cns_min_cov='${params.cns_min_cov}', run_id='${run_id}', pipeline_version='${pipeline_version}'), output_file='report.html')"
    """
    stub:
    """
    touch report.html
    """
}