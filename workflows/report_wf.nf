include { mapping_stats_table; fragment_size_table; fastp_table; kraken_table; coverage_table; rmarkdown_report } from '../modules/report'

workflow summary_report {
    take:
        consensus
        fastq_json
        kraken
        mapping_stats
        mapping_fragment_size
        mapping_coverage
        president
        pangolin
        nextclade
        nextclade_version
        nextclade_dataset_version
        sc2rf
        vois_tsv
        
    main:
        fastp_table(fastq_json.map {it -> it[1]}.collect())

        kraken_table(kraken.map {it -> it[1]}.collect(), params.taxid)

        mapping_stats_table(mapping_stats.map {it -> it[1]}.collect())

        fragment_size_table(mapping_fragment_size.map {it -> it[1]}.collect())
        
        coverage_table(mapping_coverage.map {it -> it[1]}.collect(), params.cov)

        // storeDir: "${params.output}/${params.report_dir}/"
        president_results = president.map {it -> it[1]}.collectFile(name: 'president_results.tsv', skip: 1, keepHeader: true, storeDir: "${params.output}/${params.report_dir}/single_tables")
        
        pangolin_results = pangolin.map {it -> it[1]}.collectFile(name: 'pangolin_results.tsv', skip: 1, keepHeader: true, storeDir: "${params.output}/${params.report_dir}/single_tables")

        nextclade_results = nextclade.map {it -> it[1]}.collectFile(name: 'nextclade_results.tsv', skip: 1, keepHeader: true, storeDir: "${params.output}/${params.report_dir}/single_tables")

        sc2rf_results = sc2rf.map {it -> it[1]}.collectFile(name: 'sc2rf_results.csv', skip: 1, keepHeader: true, storeDir: "${params.output}/${params.report_dir}/single_tables")
        
        vois_results = vois_tsv.map {it -> it[1]}.collectFile(name: 'vois_results.tsv', skip: 1, keepHeader: true, storeDir: "${params.output}/${params.report_dir}/single_tables")

        template = file("$baseDir/bin/summary_report.Rmd", checkIfExists: true)
        rmarkdown_report(template, fastp_table.out.stats, fastp_table.out.stats_filter, kraken_table.out.ifEmpty([]), mapping_stats_table.out, fragment_size_table.out.size, fragment_size_table.out.median, coverage_table.out.coverage_table, coverage_table.out.positive, coverage_table.out.negative, coverage_table.out.sample_cov, president_results, pangolin_results, nextclade_results, nextclade_version, nextclade_dataset_version, sc2rf_results, vois_results.ifEmpty([]))

}