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

// load adapters [optional]
if (params.adapter) { adapterInputChannel = Channel
        .fromPath( params.adapter, checkIfExists: true)
        .collect()
} else {
    adapterInputChannel = Channel
        .fromPath('NO_ADAPTERS')
        .collect()
}

// load primers [optional]
if (params.primer) { primerInputChannel = Channel
        .fromPath( params.primer, checkIfExists: true)
        .collect()
}

/************************** 
* MODULES
**************************/

// preprocess & index
include { reference_preprocessing } from './workflows/reference_preprocessing_wf'

// clip & trim
include { clip_primer } from './workflows/clip_primer_wf'
include { read_qc } from './workflows/read_qc_wf'

// read taxonomy classification
include { download_kraken_db } from './workflows/kraken_download_wf'
include { classify_reads } from './workflows/classify_reads_wf'

// map
include { mapping } from './workflows/mapping_wf'

/************************** 
* MAIN WORKFLOW
**************************/
workflow {
    // 1: reference preprocessing
    reference_preprocessing(referenceGenomeChannel)
    reference_ch = reference_preprocessing.out.ref

    // 2: quality trimming and optional adapter clipping [optional]
    reads_qc_ch = read_qc(fastqInputChannel, adapterInputChannel).reads_trimmed

    // 3: taxonomic read classification [optional]
    if (params.kraken) {
        classify_reads(reads_qc_ch, download_kraken_db())
        kraken_reports = classify_reads.out.report
    }
    reads_qc_cl_ch = params.kraken ? classify_reads.out.reads : reads_qc_ch

    // 4: read mapping
    mapping(reads_qc_cl_ch, reference_ch)

    // 5: primer clipping [optional]
    if (params.primer) {
        clip_primer(mapping.out.bam_bai, primerInputChannel)
    }

    mapping_ch = params.primer ? clip_primer.out : mapping.out.bam_bai

    // 6: variant calling
    // variants(reference_ch, reference_preprocessing.out.fai, mapping.out.bam, mapping.out.index)

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
    ${c_green}--fastq ${c_reset}                 e.g.: 'sample{1,2}.fastq' or '*.fastq.gz' or '*/*.fastq.gz'
    --list                   This flag activates csv input for the above flags [default: false]
                                 ${c_dim}style of the csv is: samplename,path_r1,path_r2${c_reset}
    --mode                          Switch between 'paired'- and 'single'-end FASTQ 

    ${c_yellow}Reference:${c_reset}
    ${c_green}--reference ${c_reset}             Currently supported: 'sars-cov2' (NC_045512)
    OR
    ${c_green}--ref_genome ${c_reset}            e.g.: 'ref.fasta'
    ${c_green}--ref_annotation ${c_reset}        e.g.: 'ref.gff'

    ${c_yellow}Adapter clipping:${c_reset}
     --adapter               Define the path of a FASTA file containing the adapter sequences to be clipped. [default: $params.adapter]

    ${c_yellow}Trimming and QC:${c_reset}
    --fastp_additional_parameters      Additional parameters for FeatureCounts [default: $params.featurecounts_additional_params]
    
    ${c_yellow}Taxonomic read filter:${c_reset}
    --kraken                 Activate taxonomic read filtering to exclude reads not classified as SARS-COV-2 (NCBI taxonomy ID 2697049) 
                                 from read mapping. A pre-processed kraken2 database will be automatically downloaded from 
                                 https://zenodo.org/record/3854856 and stored locally [default: $params.kraken]
    --taxid                  Taxonomic ID used together with the kraken2 database for read filtering [default: $params.taxid]

    ${c_yellow}Primer detection: ${c_reset}
    --primer                 Provide the path to the primer BEDPE file. [default: $params.primer]
                                 ${c_dim}TAB-delimited text file containing at least 6 fields, see here:
                                     https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format${c_reset}

    ${c_yellow}Computing options:${c_reset}
    --cores                  Max cores per process for local use [default: $params.cores]
    --max_cores              Max cores used on the machine for local use [default: $params.max_cores]
    --memory                 Max memory in GB for local use [default: $params.memory]
    --output                 Name of the result folder [default: $params.output]

    ${c_yellow}Caching:${c_reset}
    --dbs                    Location for auto-download data like databases [default: $params.dbs]
    --conda_cache_dir          Location for storing the conda environments [default: $params.conda_cache_dir]
    --singularity_cache_dir    Location for storing the singularity images [default: $params.singularity_cache_dir]
    --publish_dir_mode       Mode of output publishing: 'copy', 'symlink' [default: $params.publish_dir_mode]

    
    ${c_yellow}Execution/Engine profiles:${c_reset}
    The pipeline supports profiles to run via different ${c_green}Executers${c_reset} and ${c_blue}Engines${c_reset} e.g.: -profile ${c_green}local${c_reset},${c_blue}conda${c_reset}
    
    ${c_green}Executer${c_reset} (choose one):
      local
    
    ${c_blue}Engines${c_reset} (choose one):
      conda
    
    Per default: -profile local,conda is executed. 
    """
}