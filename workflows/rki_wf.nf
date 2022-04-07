include { rki_report } from '../modules/rki_report'

workflow rki_report_wf {
    take: 
        president_valid
        president_invalid
    main:
        readme_pdf_ch = Channel.fromPath(workflow.projectDir + "/data/rki_report/Readme.pdf")

        rki_report(president_valid.filter({ !it[0].contains("NK") }).map{it -> it[1]}.collect(), readme_pdf_ch)

        // store valid genomes
        channel_tmp1 = president_valid.filter({ !it[0].contains("NK") }).map{it -> it[1]}
        .splitText(by:100000000)
        .collectFile(name: 'valid_genomes.fasta', storeDir: params.output + "/" + params.rki_dir +"/valid/")

        // store valid genomes also as singletons
        channel_tmp2 = president_valid.filter({ !it[0].contains("NK") })
            .splitText(by:100000000)
            .collectFile(storeDir: params.output + "/" + params.rki_dir +"/valid/singletons/") { it ->
                [ "${it[0]}.fasta", it[1] ]
            }

        // store invalid genomes
        channel_tmp3 = president_invalid.filter({ !it[0].contains("NK") }).map{it -> it[1]}
            .splitText(by:100000000)
            .collectFile(name: 'invalid_genomes.fasta', storeDir: params.output + "/" + params.rki_dir +"/invalid/")

        // store invalid genomes also as singletons
        channel_tmp4 = president_invalid.filter({ !it[0].contains("NK") })
            .splitText(by:100000000)
            .collectFile(storeDir: params.output + "/" + params.rki_dir +"/invalid/singletons/") { it ->
                [ "${it[0]}.fasta", it[1] ]
            }
} 