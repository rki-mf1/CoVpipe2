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


/************************** 
* INPUT
**************************/

// load reference
if ( params.reference ) {
    if ( params.reference == 'sars-cov2' ) {
        referenceGenomeChannel = Channel
            .fromPath( workflow.projectDir + '/data/reference_SARS-CoV2/NC_045512.2.fasta' , checkIfExists: true )
            .map { file -> tuple(file.simpleName, file) }
        referenceAnnotationChannel = Channel
            .fromPath( workflow.projectDir + '/data/reference_SARS-CoV2/NC_045512.2.gff3' , checkIfExists: true )
            .map { file -> tuple(file.simpleName, file) }
    }
} else {
    referenceGenomeChannel = Channel
        .fromPath( params.ref_genome, checkIfExists: true )
        .map { file -> tuple(file.simpleName, file) }
}
if ( params.ref_annotation ) {
    referenceAnnotationChannel = Channel
        .fromPath( params.ref_annotation, checkIfExists: true )
        .map { file -> tuple(file.simpleName, file) }
}

// illumina reads input & --list support
if (params.mode == 'paired') {
    if (params.fastq && params.list) { fastqInputChannel = Channel
        .fromPath( params.fastq, checkIfExists: true )
        .splitCsv()
        .map { row -> [row[0], [file(row[1], checkIfExists: true), file(row[2], checkIfExists: true)]] }
        }
    else if (params.fastq) { fastqInputChannel = Channel
        .fromFilePairs( params.fastq, checkIfExists: true )
        }
} else {
    if (params.fastq && params.list) { fastqInputChannel = Channel
        .fromPath( params.fastq, checkIfExists: true )
        .splitCsv()
        .map { row -> [row[0], [file(row[1], checkIfExists: true)]] }
        }
    else if (params.fastq) { fastqInputChannel = Channel
        .fromPath( params.fastq, checkIfExists: true )
        .map { file -> [file.simpleName, [file]]}
        }
}

/************************** 
* MODULES
**************************/

// preprocess & index
include { dos2unix } from './modules/dos2unix'
include { index_samtools } from './modules/samtools'
include { index_picard } from './modules/picard'
include { index_bwa } from './modules/bwa'


/************************** 
* DATABASES
**************************/

/************************** 
* SUB WORKFLOWS
**************************/

// 1: indexing
workflow indexing {
    take: reference_fasta
    main:
        dos2unix(reference_fasta)

        index_samtools(dos2unix.out)
        index_picard(dos2unix.out) // this seems to be unused
        index_bwa(dos2unix.out)

    emit: 
        ref = dos2unix.out
        fai = index_samtools.out
        bwa = index_bwa.out.index

}

/************************** 
* MAIN WORKFLOW
**************************/
workflow {

    // generate all indices
    indexing(referenceGenomeChannel)
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