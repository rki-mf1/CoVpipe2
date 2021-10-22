#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// help message
if (params.help) { exit 0, helpMSG() }

// error codes
if (params.profile) {
    exit 1, "--profile is WRONG use -profile" }
if (params.workdir) {
    exit 1, "--workdir is WRONG use -w" }

// warnings
if ( workflow.profile == 'standard' ) { 
    "NO EXECUTION PROFILE SELECTED, using [-profile local,conda]" }

if ( !workflow.revision ) { 
    println ""
    println "\u001B[31mWARN: no revision selected, please use -r for full reproducibility\033[0m"
}


/************************** 
* PARAMETERS
**************************/

if ( !params.fastq ) {
    exit 1, "input missing, use [--fastq]"
}

Set reference = ['sars-cov2'] // can be extended later on
if ( !params.reference && !params.ref_genome && !params.ref_annotation ) {
    exit 1, "reference missing, use [--ref_genome] (and [--ref_annotation]) or choose of " + reference + " with [--reference]"
}
if ( params.reference && params.ref_genome ) {
    exit 1, "too many references, use either [--ref_genome] (and [--ref_annotation]), or [--reference]"
}
if ( params.reference && ! (params.reference in reference) ) {
    exit 1, "unknown reference, currently supported: " + reference
}
// kraken input test
if (params.kraken && ! params.taxid) {
    exit 1, "Kraken2 database defined but no --taxid!"
}

/************************** 
* INPUT
**************************/

// load reference
if ( params.reference ) {
    if ( params.reference == 'sars-cov2' ) {
        referenceGenomeChannel = Channel
            .fromPath( workflow.projectDir + '/data/reference_SARS-CoV2/NC_045512.2.fasta' , checkIfExists: true )
        referenceAnnotationChannel = Channel
            .fromPath( workflow.projectDir + '/data/reference_SARS-CoV2/NC_045512.2.gff3' , checkIfExists: true )
    }
} else {
    referenceGenomeChannel = Channel
        .fromPath( params.ref_genome, checkIfExists: true )
}
if ( params.ref_annotation ) {
    referenceAnnotationChannel = Channel
        .fromPath( params.ref_annotation, checkIfExists: true )
}

// illumina reads input & --list support
if (params.mode == 'paired') {
    if (params.fastq && params.list) { fastqInputChannel = Channel
        .fromPath( params.fastq, checkIfExists: true )
        .splitCsv()
        .map { row -> [row[0], [file(row[1], checkIfExists: true), file(row[2], checkIfExists: true)]] }}
    else if (params.fastq) { fastqInputChannel = Channel
        .fromFilePairs( params.fastq, checkIfExists: true )}
} else {
    if (params.fastq && params.list) { fastqInputChannel = Channel
        .fromPath( params.fastq, checkIfExists: true )
        .splitCsv()
        .map { row -> [row[0], [file(row[1], checkIfExists: true)]] }}
    else if (params.fastq) { fastqInputChannel = Channel
        .fromPath( params.fastq, checkIfExists: true )
        .map { file -> [file.simpleName, [file]]}}
}

// load primers [optional]
if (params.primer) { primerInputChannel = Channel
        .fromPath( params.primer, checkIfExists: true)
        .collect()
}

// load adapters [optional]
if (params.adapter) { adapterInputChannel = Channel
        .fromPath( params.adapter, checkIfExists: true)
        .collect()
} else {
    adapterInputChannel = Channel
        .fromPath('NO_ADAPTERS')
        .collect()
}

/************************** 
* MODULES
**************************/

// preprocess & index
include { dos2unix } from './modules/dos2unix'

// clip & trim
include { trim_primer } from './modules/trimprimer'
include { fastp } from './modules/fastp'

// read taxonomy classification
include { kraken_db; kraken; filter_virus_reads } from './modules/kraken'

// map
include { index_bwa; bwa } from './modules/bwa'
include { get_genomecov } from './modules/bedtools'

// utils
include { index_fasta; index_bam } from './modules/samtools'

/************************** 
* DATABASES
**************************/

/* Comment section:
The Database Section is designed to "auto-get" pre-prepared databases.
It is written for local use and perspective cloud use via params.cloudProcess (see nextflow.config).
*/

workflow download_kraken_db {
    main:
        if (params.kraken) {
            // local storage via storeDir
            if (!params.cloudProcess) { kraken_db() ; database_kraken = kraken_db.out}
            // cloud storage file.exists()?
            if (params.cloudProcess) { 
                kraken_db_preload = file("${params.databases}/kraken/GRCh38.p13_GBcovid19-2020-05-22.tar.gz")
                if (kraken_db_preload.exists()) { database_kraken = kraken_db_preload }    
                else { kraken_db() ; database_kraken = kraken_db.out }
            }
        }
    emit: database_kraken.collect()
}  


/************************** 
* SUB WORKFLOWS
**************************/

// 1: reference preprocessing
workflow reference_preprocessing {
    take: reference_fasta
    main:
        dos2unix(reference_fasta) \
            | index_fasta
    emit: 
        ref = dos2unix.out
        fai = index_fasta.out
}

// 2: amplicon primer clipping [optional]
workflow primer {
    take:
        illumina_reads
        primer_set
    main:
        trim_primer(illumina_reads, primer_set)
    emit:
        trim_primer.out.reads
}

// 3: quality trimming and optional adapter clipping [optional]
workflow read_qc {
    take: 
        illumina_reads
        adapter_fasta
    main:
        fastp(illumina_reads, adapter_fasta)
        reads_trimmed = fastp.out.reads
        fastp_json = fastp.out.json
    emit: 
        reads_trimmed
        fastp_json
}

// 4: taxonomic read classification [optional]
workflow classify {
    take:
        reads
        db
    main:
        filter_virus_reads(kraken(reads, db).fastq)
    emit:
        reads = filter_virus_reads.out.fastq
        report = kraken.out.kraken_report
}

// 5: read mapping
workflow mapping {
    take: 
        illumina_reads
        reference_fasta
    main:

        index_bwa(reference_fasta)
        bwa(illumina_reads, index_bwa.out.collect()) |
            (index_bam & get_genomecov)
    emit:
        bam = bwa.out.bam
        index = index_bam.out.index
        coverage = get_genomecov.out.tsv
}


/************************** 
* MAIN WORKFLOW
**************************/
workflow {

    // generate all indices for the reference
    reference_preprocessing(referenceGenomeChannel)
    reference_ch = reference_preprocessing.out.ref

    // primer clipping [optional]
    if (params.primer) {
        primer(fastqInputChannel, primerInputChannel)
    }
    reads_ch = params.primer ? primer.out : fastqInputChannel

    // quality and adapter trimming
    reads_qc_ch = read_qc(reads_ch, adapterInputChannel).reads_trimmed

    // taxonomic read classification
    if (params.kraken) {
        classify(reads_qc_ch, download_kraken_db())
        kraken_reports = classify.out.report
    }
    reads_qc_cl_ch = params.kraken ? classify.out.reads : reads_qc_ch

    // read mapping
    mapping(reads_qc_cl_ch, reference_ch)

}

/************************** 
* HELP
**************************/
def helpMSG() {
    c_green = "\033[0;32m";
    c_reset = "\033[0m";
    c_yellow = "\033[0;33m";
    c_blue = "\033[0;34m";
    c_red = "\u001B[31m";
    c_dim = "\033[2m";
    log.info """
    ____________________________________________________________________________________________
    
    ${c_blue}Robert Koch Institute, MF1 Bioinformatics${c_reset}

    Workflow: CoVpipeNext

    ${c_yellow}Usage examples:${c_reset}
    nextflow run CoVpipeNext.nf --fastq '*R{1,2}.fastq.gz' --reference 'sars-cov2' --cores 4 --max_cores 8
    or
    nextflow run RKIBioinformaticsPipelines/covpipenxt -r <version> --fastq '*R{1,2}.fastq.gz' --reference ref.fasta --cores 4 --max_cores 8

    ${c_yellow}Inputs:
    Illumina read data:${c_reset}
    ${c_green}--fastq ${c_reset}            e.g.: 'sample{1,2}.fastq' or '*.fastq.gz' or '*/*.fastq.gz'
    --list              this flag activates csv input for the above flags [default: false]
                        style of the csv is: ${c_dim}samplename,path_r1,path_r2${c_reset}
    --mode                     switch between 'paired'- and 'single'-end FASTQ 
    ${c_yellow}Reference:${c_reset}
    ${c_green}--reference ${c_reset}        currently supported: SARS-CoV2 (NC_045512)
    OR
    ${c_green}--ref_genome ${c_reset}       e.g.: 'ref.fasta'
    ${c_green}--ref_annotation ${c_reset}   e.g.: 'ref.fasta'
    """
}