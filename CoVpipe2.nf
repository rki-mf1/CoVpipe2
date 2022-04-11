#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// help message
if (params.help) { exit 0, helpMSG() }

// parameter sanity check
Set valid_params = ['cores', 'max_cores', 'memory', 'help', 'profile', 'workdir', 'fastq', 'list', 'dir', 'mode', 'run_id', 'reference', 'ref_genome', 'ref_annotation', 'adapter', 'fastp_additional_parameters', 'kraken', 'taxid', 'primer_bed', 'primer_bedpe', 'primer_version', 'vcount', 'frac', 'cov', 'vois', 'var_mqm', 'var_sap', 'var_qual', 'cns_min_cov', 'cns_gt_adjust', 'update', 'pangolin_docker_default', 'nextclade_docker_default', 'output', 'reference_dir', 'read_dir', 'mapping_dir', 'variant_calling_dir', 'consensus_dir', 'linage_dir', 'report_dir', 'rki_dir', 'runinfo_dir', 'singularity_cache_dir', 'conda_cache_dir', 'databases', 'publish_dir_mode', 'cloudProcess', 'cloud-process']
def parameter_diff = params.keySet() - valid_params
if (parameter_diff.size() != 0){
    exit 1, "ERROR: Parameter(s) $parameter_diff is/are not valid in the pipeline!\n"
}

// error codes
if (params.profile) {
    exit 1, "--profile is WRONG use -profile" }
if (params.workdir) {
    exit 1, "--workdir is WRONG use -w" }

// warnings
def folder = new File(params.output)
if ( folder.exists() ) { 
    println ""
    println "\033[0;33mWARNING: Output folder already exists. Results might be overwritten! You can adjust the output folder via [--output]\033[0m"
}
if ( !workflow.revision ) { 
    println ""
    println "\033[0;33mWARNING: not a stable execution. Please use -r for full reproducibility.\033[0m"
}
if ( workflow.profile.contains('singularity') ) {
    println ""
    println "\033[0;33mWARNING: Singularity image building sometimes fails!"
    println "Multiple resumes (-resume) and --max_cores 1 --cores 1 for local execution might help.\033[0m\n"
}

// print info message
defaultMSG()

/************************** 
* PARAMETERS
**************************/

if ( !params.fastq ) {
    exit 1, "input missing, use [--fastq]"
}

Set reference = ['sars-cov-2'] // can be extended later on
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
    if ( params.reference == 'sars-cov-2' ) {
        ref_genome_file = file( workflow.projectDir + '/data/reference/SARS-CoV-2/MN908947.3.fasta' , checkIfExists: true )
        ref_annotation_file = file( workflow.projectDir + '/data/reference/SARS-CoV-2/MN908947.3.gff3' , checkIfExists: true )
    }
} else {
    if ( params.ref_genome ) {
        ref_genome_file = file( params.ref_genome, checkIfExists: true )
    }
    if ( params.ref_annotation ) {
        ref_annotation_file = file( params.ref_annotation, checkIfExists: true )
    }
}

// illumina reads input & --list support
if (! params.dir) {
    if (params.mode == 'paired') {
        if (params.fastq && params.list) { fastqInputChannel = Channel
            .fromPath( params.fastq, checkIfExists: true )
            .splitCsv(header: true, sep: ',')
            .map { row -> [row.sample, [file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true)]] }
        } else if (params.fastq) { fastqInputChannel = Channel
            .fromFilePairs( params.fastq, checkIfExists: true )}
    } else {
        if (params.fastq && params.list) { fastqInputChannel = Channel
            .fromPath( params.fastq, checkIfExists: true )
            .splitCsv(header: true, sep: ',')
            .map { row -> [row.sample, [file(row.fastq, checkIfExists: true)]] }}
        else if (params.fastq) { fastqInputChannel = Channel
            .fromPath( params.fastq, checkIfExists: true )
            .map { file -> [file.simpleName, [file]]}}
    }
}

// load adapters [optional]
adapter_file = params.adapter ? file(params.adapter, checkIfExists: true) : file('NO_ADAPTERS')

// load primers [optional]
list = [params.primer_bedpe, params.primer_bed, params.primer_version]
assert list.count(false) >= 2: "Choose either one of these tree parameters (--primer_bedpe, --primer_bed, --primer_version) or none of them."
if( params.primer_bedpe || params.primer_bed || params.primer_version ){

    if( params.primer_bedpe ){
        primer_file = file(params.primer_bedpe, checkIfExists: true)
    } else if ( params.primer_bed ){
        primer_file = file(params.primer_bed, checkIfExists: true)
    } else if ( params.primer_version ){
        primer_file = file(workflow.projectDir + "/data/external_primer_schemes/SARS-CoV-2/${params.primer_version}/nCoV-2019.scheme.bed", checkIfExists: true)
    } else {
        println "No primer input."
    }

    // check if fasta header matches with primer chrom
    // if not exit, because it does not do what the user expects
    // will fail for multi fastas!
    primer_file.withReader { line = it.readLine() } // read only first line
    primer_id = line.split()[0] // extract id

    ref_genome_file.withReader { line = it.readLine() }  // read only first line
    ref_id = line.split()[0].replaceAll("^>", "") // extract id

    assert ref_id == primer_id: "Faster header ($ref_id) and primer chrom ($primer_id) don't match. Provide a matching primer BEDPE file or don't set --primer"

    primer_file_ch = Channel.fromPath(primer_file, checkIfExists: true)
}

// load vois [optional]
if( params.vois ){ 
    vois_file = file(params.vois, checkIfExists: true) 

    primer_file.withReader { line = it.readLine() } // read only first line
    primer_id = line.split()[0] // extract id

    voi_id=vois_file.withReader {  r-> r.eachLine {it} }.split()[0]
    assert ref_id == voi_id: "Faster header ($ref_id) and VOI VCF file chrom ($voi_id) don't match. Provide a matching VCF file or don't set --vois"
}

/************************** 
* AUTO-UPDATE (nextclade and pangolin)
**************************/

static boolean DockernetIsAvailable() {
    try {
        final URL url = new URL("https://registry.hub.docker.com/v2/repositories/rkimf1/pangolin/tags/");
        final URLConnection conn = url.openConnection();
        conn.connect();
        conn.getInputStream().close();
        return true;
    } catch (MalformedURLException e) {
        return false;
    } catch (IOException e) {
        return false;
    }
}

def internetcheck = DockernetIsAvailable()

if ( ( workflow.profile.contains('singularity') || workflow.profile.contains('docker') ) && params.update) {
println "\033[0;33mWarning: Running --update might not be CoVpipe compatible!\033[0m"
    if ( internetcheck.toString() == "true" ) { 
        tagname = 'https://registry.hub.docker.com/v2/repositories/rkimf1/pangolin/tags/'.toURL().text.split(',"name":"')[1].split('","')[0]
        params.pangolin_docker = "rkimf1/pangolin:" + tagname
        println "\033[0;32mFound latest pangolin container, using: " + params.pangolin_docker + " \033[0m" 

        tagname = 'https://registry.hub.docker.com/v2/repositories/rkimf1/nextclade/tags/'.toURL().text.split(',"name":"')[1].split('","')[0]
        params.nextclade_docker = "rkimf1/nextclade:" + tagname 
        println "\033[0;32mFound latest nextclade container, using: " + params.nextclade_docker + " \033[0m"
    } 
    if ( internetcheck.toString() == "false" ) { 
        println "\033[0;33mCould not find the latest pangolin container, trying: " + params.pangolin_docker_default + "\033[0m"
        params.pangolin_docker = params.pangolin_docker_default 

        println "\033[0;33mCould not find the latest nextclade container, trying: " + params.nextclade_docker_default + "\033[0m"
        params.nextclade_docker = params.nextclade_docker_default 
    } 
}
else { params.pangolin_docker = params.pangolin_docker_default ; params.nextclade_docker = params.nextclade_docker_default  }

/************************** 
* MODULES
**************************/

// preprocess & index
include { reference_preprocessing } from './workflows/reference_preprocessing_wf'

// qc & trim
include { read_qc } from './workflows/read_qc_wf'

// read taxonomy classification
include { download_kraken_db } from './workflows/kraken_download_wf'
include { classify_reads } from './workflows/classify_reads_wf'

// map
include { mapping } from './workflows/mapping_wf'

// primer clipping
include { clip_primer } from './workflows/clip_primer_wf'

// variant calling and annotation
include { variant_calling } from './workflows/variant_calling_wf'
include { annotate_variant } from './workflows/annotate_variant_wf'

// generate consensus
include { generate_consensus } from './workflows/generate_consensus_wf'
include { annotate_consensus } from './workflows/annotate_consensus_wf'

// variants of interest
include { inspect_vois } from './workflows/inspect_vois_wf'

// assign linages
include { assign_linages } from './workflows/assign_linages_wf'
// genome quality (president)
include { genome_quality } from './workflows/genome_quality_wf'

include { summary_report } from './workflows/report_wf'
include { rki_report_wf } from './workflows/rki_wf'

include { bed2bedpe; make_sample_sheet } from './modules/utils'

/************************** 
* MAIN WORKFLOW
**************************/
workflow {
    if(params.dir){
        make_sample_sheet(params.fastq)

        if (params.mode == 'paired') { fastqInputChannel = make_sample_sheet.out
            .splitCsv(header: true, sep: ',')
            .map { row -> [row.sample, [file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true)]] }
        } else { fastqInputChannel = make_sample_sheet.out
            .splitCsv(header: true, sep: ',')
            .map { row -> [row.sample, [file(row.fastq, checkIfExists: true)]] }
        }
    }

    // 1: reference preprocessing
    reference_preprocessing(ref_genome_file)
    reference_ch = reference_preprocessing.out.ref

    // 2: quality trimming and optional adapter clipping
    reads_qc_ch = read_qc(fastqInputChannel, adapter_file).reads_trimmed

    // 3: taxonomic read classification [optional]
    if (params.kraken) {
        classify_reads(reads_qc_ch, download_kraken_db())
        kraken_reports = classify_reads.out.report
    } else {
        kraken_reports = Channel.empty()
    }
    reads_qc_cl_ch = params.kraken ? classify_reads.out.reads : reads_qc_ch

    // 4: read mapping
    mapping(reads_qc_cl_ch, reference_ch)

    // 5: primer clipping [optional]
    if (params.primer_version || params.primer_bedpe || params.primer_bed) {
        if( params.primer_version || params.primer_bed ){
            new_basename = params.primer_version ? params.primer_version : primer_file.baseName
            bed2bedpe(primer_file_ch.map{primer_scheme -> [new_basename, primer_scheme]}, '_LEFT', '_RIGHT')
            primer_file_ch = bed2bedpe.out
        }
        clip_primer(mapping.out.bam_bai, primer_file_ch)
    }
    mapping_ch = params.primer_version || params.primer_bedpe || params.primer_bed ? clip_primer.out : mapping.out.bam_bai

    // 6: variant calling
    variant_calling(reference_ch, reference_preprocessing.out.fai, mapping_ch)

    // 7: generate consensus
    generate_consensus(variant_calling.out.vcf, reference_ch, mapping_ch.map{ it[0,1]})

    // 8: annotate consensus [optional]
    if (params.reference || params.ref_annotation) {
        annotate_consensus(generate_consensus.out.consensus_ambiguous.mix(generate_consensus.out.consensus_masked), reference_ch, ref_annotation_file)
    }

    // 9: annotate mutations
    annotate_variant(variant_calling.out.vcf, generate_consensus.out.consensus_ambiguous, 'sars-cov-2', 'NC_045512.2')

    // 10: compare with variants/mutations of interest [optional]
    if ( params.vois ) {
        inspect_vois(vois_file, variant_calling.out.vcf_csi, generate_consensus.out.low_coverage_bed)
        vois = inspect_vois.out
    } else {
        vois = Channel.empty()
    }

    // 11: linage assignment, genome quality
    assign_linages(generate_consensus.out.consensus_ambiguous)
    genome_quality(generate_consensus.out.consensus_ambiguous, reference_ch)

    // 12: report
    summary_report(generate_consensus.out.consensus_ambiguous, read_qc.out.fastp_json, kraken_reports.ifEmpty([]), mapping.out.mapping_stats, mapping.out.fragment_size, mapping.out.coverage, genome_quality.out.report, assign_linages.out.report, annotate_variant.out.nextclade_results, annotate_variant.out.nextclade_version, annotate_variant.out.nextclade_dataset_version, vois.ifEmpty([]) )

    // 13: provide data for DESH upload at RKI
    rki_report_wf(genome_quality.out.valid, genome_quality.out.invalid)
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

    Workflow: CoVpipe2

    ${c_yellow}Usage examples:${c_reset}
    nextflow run CoVpipe2.nf --fastq '*R{1,2}.fastq.gz' --reference 'sars-cov-2' --cores 4 --max_cores 8
    or
    nextflow run rki-mf1/CoVpipe2 -r <version> --fastq '*R{1,2}.fastq.gz' --ref_genome ref.fasta --cores 4 --max_cores 8

    ${c_yellow}Reference, required:${c_reset}
    ${c_green}--reference ${c_reset}             Currently supported: 'sars-cov-2' (MN908947.3)
    OR
    ${c_green}--ref_genome ${c_reset}            Reference FASTA file.
    ${c_green}--ref_annotation ${c_reset}        Reference GFF file.

    ${c_yellow}Illumina read data, required:${c_reset}
    ${c_green}--fastq ${c_reset}                 e.g.: 'sample{1,2}.fastq' or '*.fastq.gz' or '*/*.fastq.gz'

    ${c_yellow}Optional input settings:${c_reset}
    --list                   This flag activates csv input for --fastq [default: false]
                                 ${c_dim}style and header of the csv is: sample,fastq_1,fastq_2${c_reset}
    --mode                   Switch between 'paired'- and 'single'-end FASTQ; 'single' is experimental [default: $params.mode]
    --run_id                 Run ID [default: $params.run_id]

    ${c_yellow}Adapter clipping:${c_reset}
     --adapter               Define the path of a FASTA file containing the adapter sequences to be clipped. [default: $params.adapter]

    ${c_yellow}Trimming and QC:${c_reset}
    --fastp_additional_parameters      Additional parameters for FeatureCounts [default: $params.fastp_additional_parameters]
    
    ${c_yellow}Taxonomic read filter:${c_reset}
    --kraken                 Activate taxonomic read filtering to exclude reads not classified with specific taxonomic ID (see --taxid) [default: $params.kraken]
                                 ${c_dim}A pre-processed kraken2 database will be automatically downloaded from 
                                 https://zenodo.org/record/3854856 and stored locally.${c_reset}
    --taxid                  Taxonomic ID used together with the kraken2 database for read filtering [default: $params.taxid]

    ${c_yellow}Primer detection: ${c_reset}
    --primer_bedpe           Provide the path to the primer BEDPE file. [default: $params.primer_bedpe]
                                 ${c_dim}TAB-delimited text file containing at least 6 fields, see here:
                                 https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format${c_reset}
    OR
    --primer_bed             Provide the path to the primer BED file. [default: $params.primer_bed]
    OR
    --primer_version         Provide a primer version. Currently supported ARTIC versions: V1, V2, V3, V4, V4.1 [default: $params.primer_version]

    ${c_yellow}Variant calling:${c_reset}
    --vcount                 Minimum number of reads at a position to be considered for variant calling. [default: $params.vcount]
    --cov                    Minimum number of supporting reads which are required to call a variant. [default: $params.cov]
    --frac                   Minimum percentage of supporting reads at the respective position required to call a variant. 
                                 In turn, variants supported by (1 - frac)*100% reads will be explicitly called. [default: $params.frac]
    --vois                   Compare called variants to a VCF file with you variants of interest [default: $params.vois]

    ${c_yellow}Variant hard filtering:${c_reset}
    --var_mqm                Minimal mean mapping quality of observed alternate alleles (MQM). The mapping quality (MQ) 
                                 measures how good reads align to the respective reference genome region. Good mapping qualities are 
                                 around MQ 60. GATK recommends hard filtering of variants with MQ less than 40. [default: $params.var_mqm]
    --var_sap                Strand balance probability for the alternate allele (SAP). The SAP is the Phred-scaled 
                                 probability that there is strand bias at the respective site. A value near 0 indicates little or 
                                 no strand bias. Set to -1 to disable the filter. [default: $params.var_sap]
    --var_qual               Minimal variant call quality. Freebayes produces a general judgement of the 
                                 variant call. [default: $params.var_qual]

    ${c_yellow}Consensus generation:${c_reset}
    --cns_min_cov            Minimum number of reads required so that the respective position in the consensus sequence 
                                 is NOT hard masked. [default: $params.cns_min_cov]
    --cns_gt_adjust          Minimum fraction of reads supporting a variant which leads to an explicit call of this 
                                 variant (genotype adjustment). The value has to be greater than 0.5 but not greater than 1. 
                                 To turn genotype adjustment off, set the value to 0. [default: $params.cns_gt_adjust]

    ${c_yellow}Updated for linage assignment and mutation calling:${c_reset}
    --update                   Update pangolin and nextclade [default: $params.update]
                                  ${c_dim}Depending on the chosen engine either the conda environment (profiles: 'standard', 'conda', 'mamba') 
                                  or the container (profiles: 'docker', 'singularity') is updated.${c_reset}
    --pangolin_docker_default  Default container tag for pangolin [default: $params.pangolin_docker_default]
    --nextclade_docker_default Default container tag for nextclade [default: $params.nextclade_docker_default]

    ${c_yellow}Computing options:${c_reset}
    --cores                  Max cores per process for local use [default: $params.cores]
    --max_cores              Max cores used on the machine for local use [default: $params.max_cores]
    --memory                 Max memory in GB for local use [default: $params.memory]

    ${c_yellow}Output options:${c_reset}
    --output                 Name of the result folder [default: $params.output]
    --publish_dir_mode       Mode of output publishing: 'copy', 'symlink' [default: $params.publish_dir_mode]
                                 ${c_dim}With 'symlink' results are lost when removing the work directory.${c_reset}

    ${c_yellow}Caching:${c_reset}
    --databases              Location for auto-download data like databases [default: $params.databases]
    --conda_cache_dir        Location for storing the conda environments [default: $params.conda_cache_dir]
    --singularity_cache_dir  Location for storing the singularity images [default: $params.singularity_cache_dir]
    
    ${c_yellow}Execution/Engine profiles:${c_reset}
    The pipeline supports profiles to run via different ${c_green}Executers${c_reset} and ${c_blue}Engines${c_reset} e.g.: -profile ${c_green}local${c_reset},${c_blue}conda${c_reset}
    
    ${c_green}Executer${c_reset} (choose one):
      local
      slurm
    
    ${c_blue}Engines${c_reset} (choose one):
      conda
      mamba
      docker
      singularity

    ${c_dim}Misc:
      cluster                Loads resource configs more suitable for cluster execution.
                             Has to be combine with an engine and an executor.
    ${c_reset}

    Per default: -profile local,conda is executed. 
    """
}
def defaultMSG(){
    println " "
    println "\u001B[32mProfile: $workflow.profile\033[0m"
    println " "
    println "\033[2mCurrent User: $workflow.userName"
    println "Nextflow-version: $nextflow.version"
    println "Starting time: $workflow.start"
    println "Workdir location:"
    println "  $workflow.workDir"
    println "Launchdir location:"
    println "  $workflow.launchDir"
    println "Permanent cache directory:"
    println "  $params.databases"
    if ( workflow.profile.contains('conda') || workflow.profile.contains('mamba') || workflow.profile.contains('standard') ) { 
        println "Conda/mamba cache directory:"
        println "  $params.conda_cache_dir"
    }
    println "Configuration files:"
    println "  $workflow.configFiles"
    println "Cmd line:"
    println "  $workflow.commandLine\u001B[0m"
    if (workflow.repository != null){ println "\033[2mGit info: $workflow.repository - $workflow.revision [$workflow.commitId]\u001B[0m" }
    println " "
    if (workflow.profile.contains('standard') || workflow.profile.contains('local')) {
        println "\033[2mCPUs to use: $params.cores, maximal CPUs to use: $params.max_cores\u001B[0m"
    }
}