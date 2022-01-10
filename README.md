# CoVpipeNext

## Quick installation

The pipeline is written in [`Nextflow`](https://nf-co.re/usage/installation), which can be used on any POSIX compatible system (Linux, OS X, etc). Windows system is supported through [WSL](https://en.wikipedia.org/wiki/Windows_Subsystem_for_Linux). You need `Nextflow` installed and either `conda` to run the steps of the pipeline:

1. Install  `Nextflow`
    <details><summary>click here for a bash one-liner </summary>

    ```bash
    wget -qO- https://get.nextflow.io | bash
    # In the case you don’t have wget
    # curl -s https://get.nextflow.io | bash
    ```

    </details>
2. Install [`conda`](https://conda.io/miniconda.html)
    <details><summary>click here for a bash two-liner for Miniconda3 Linux 64-bit</summary>

    ```bash
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh
    ```

    </details>

OR

1. Install `conda`
    <details><summary>click here for a bash two-liner for Miniconda3 Linux 64-bit</summary>

    ```bash
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh
    ```

    </details>
1. Install `Nextflow` via `conda`
    <details><summary>click here to see how to do that</summary>

    ```bash
    conda create -n nextflow -c bioconda nextflow
    conda active nextflow
    ```

    </details>

All other dependencies and tools will be installed within the pipeline via `conda`.

### Call help

```bash
nextflow run RKIBioinformaticsPipelines/covpipenext -hub gitlab --help
```

### Update the pipeline

```bash
nextflow pull RKIBioinformaticsPipelines/covpipenext -hub gitlab
```

### Use a certain release

We recommend to use a of the pipeline:

```bash
nextflow pull RKIBioinformaticsPipelines/covpipenext -hub gitlab -r <RELEASE>
```

## Help message

<details><summary>click here to see the complete help message</summary>

```
    Robert Koch Institute, MF1 Bioinformatics

    Workflow: CoVpipeNext

    Usage examples:
    nextflow run CoVpipeNext.nf --fastq '*R{1,2}.fastq.gz' --reference 'sars-cov2' --cores 4 --max_cores 8
    or
    nextflow run RKIBioinformaticsPipelines/covpipenxt -r <version> --fastq '*R{1,2}.fastq.gz' --reference ref.fasta --cores 4 --max_cores 8

    Inputs:
    Illumina read data:
    --fastq                  e.g.: 'sample{1,2}.fastq' or '*.fastq.gz' or '*/*.fastq.gz'
    --list                   This flag activates csv input for the above flags [default: false]
                                 style and header of the csv is: samplename,path_r1,path_r2
    --mode                          Switch between 'paired'- and 'single'-end FASTQ [default: paired]
    --run_id                 Run ID [default: ]

    Reference:
    --reference              Currently supported: 'sars-cov2' (NC_045512)
    OR
    --ref_genome             e.g.: 'ref.fasta'
    --ref_annotation         e.g.: 'ref.gff'

    Adapter clipping:
     --adapter               Define the path of a FASTA file containing the adapter sequences to be clipped. [default: false]

    Trimming and QC:
    --fastp_additional_parameters      Additional parameters for FeatureCounts [default: --qualified_quality_phred 20 --length_required 50]
    
    Taxonomic read filter:
    --kraken                 Activate taxonomic read filtering to exclude reads not classified as SARS-COV-2 (NCBI taxonomy ID 2697049) 
                                 from read mapping. A pre-processed kraken2 database will be automatically downloaded from 
                                 https://zenodo.org/record/3854856 and stored locally [default: false]
    --taxid                  Taxonomic ID used together with the kraken2 database for read filtering [default: 2697049]

    Primer detection: 
    --primer                 Provide the path to the primer BEDPE file. [default: false]
                                 TAB-delimited text file containing at least 6 fields, see here:
                                 https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format

    Variant calling:
    --vcount                 Minimum number of reads at a position to be considered for variant calling. [default: 10]
    --cov                    Minimum number of supporting reads which are required to call a variant. [default: 20]
    --frac                   Minimum percentage of supporting reads at the respective position required to call a variant. 
                                 In turn, variants supported by (1 - frac)*100% reads will be explicitly called. [default: 0.1]
    --vois                   Compare called variants to a VCF file with you variants of interest [default: false]

    Variant hard filtering:
    --var_mqm                Minimal mean mapping quality of observed alternate alleles (MQM). The mapping quality (MQ) 
                                measures how good reads align to the respective reference genome region. Good mapping qualities are 
                                around MQ 60. GATK recommends hard filtering of variants with MQ less than 40. [default: 40]
    --var_sap                Strand balance probability for the alternate allele (SAP). The SAP is the Phred-scaled 
                                probability that there is strand bias at the respective site. A value near 0 indicates little or 
                                no strand bias.  [default: 60]
    --var_qual               Minimal variant call quality. Freebayes produces a general judgement of the 
                                variant call. [default: 10]

    Consensus generation:
    --cns_min_cov            Minimum number of reads required so that the respective position in the consensus sequence 
                                 is NOT hard masked. [default: 20]
    --cns_gt_adjust          Minimum fraction of reads supporting a variant which leads to an explicit call of this 
                                 variant (genotype adjustment). The value has to be greater than 0.5 but not greater than 1. 
                                 To turn genotype adjustment off, set the value to 0. [default: 0.9]

    Linage assignment:
    --update_pangolin        Update pangolin environment to get the latest version that is available from bioconda.

    Computing options:
    --cores                  Max cores per process for local use [default: 4]
    --max_cores              Max cores used on the machine for local use [default: 12]
    --memory                 Max memory in GB for local use [default: 12]
    --output                 Name of the result folder [default: results]

    Caching:
    --databases                Location for auto-download data like databases [default: nextflow-autodownload-databases]
    --conda_cache_dir          Location for storing the conda environments [default: conda]
    --singularity_cache_dir    Location for storing the singularity images [default: singularity]
    --publish_dir_mode       Mode of output publishing: 'copy', 'symlink' [default: copy]

    
    Execution/Engine profiles:
    The pipeline supports profiles to run via different Executers and Engines e.g.: -profile local,conda
    
    Executer (choose one):
      local
      slurm
    
    Engines (choose one):
      conda

    Misc:
      cluster                Loads resource configs more suitable for cluster execution.
                             Has to be combine with an engine and an executor.
    

    Per default: -profile local,conda is executed. 
```

</details>