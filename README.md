# CoVpipe2

CoVpipe2 is a Nextflow pipeline for reference-based genome reconstruction of SARS-CoV-2 from NGS data. In principle it can be used also for other viruses.

<details><summary>Table of contents</summary>

- [CoVpipe2](#covpipe2)
  - [Quick installation](#quick-installation)
    - [Call help](#call-help)
    - [Update the pipeline](#update-the-pipeline)
    - [Use a certain release](#use-a-certain-release)
  - [Quick run examples](#quick-run-examples)
    - [Example 1:](#example-1)
    - [Example 2:](#example-2)
    - [Example sample sheet](#example-sample-sheet)
  - [Manual / help](#manual--help)
  - [Workflow](#workflow)
  - [Acknowledgement, props and inspiration](#acknowledgement-props-and-inspiration)

</details>

## Quick installation

The pipeline is written in [`Nextflow`](https://nf-co.re/usage/installation), which can be used on any POSIX compatible system (Linux, OS X, etc). Windows system is supported through [WSL](https://en.wikipedia.org/wiki/Windows_Subsystem_for_Linux). You need `Nextflow` installed and either `conda`, or `Docker`, or `Singularity` to run the steps of the pipeline:

- Install  `Nextflow` via self-installing package
    <details><summary>click here for a bash one-liner </summary>

    ```bash
    wget -qO- https://get.nextflow.io | bash
    # In the case you don’t have wget
    # curl -s https://get.nextflow.io | bash
    ```

    </details>

OR

- Install `Nextflow` via `conda`
    <details><summary>click here for a bash two-liner for Miniconda3 Linux 64-bit</summary>

    ```bash
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh
    conda create -n nextflow -c bioconda nextflow
    conda active nextflow
    ```

    </details>

All other dependencies and tools will be installed within the pipeline via `conda`, `Docker` or `Singularity`.

### Call help

```bash
nextflow run rki-mf1/CoVpipe2 --help
```

### Update the pipeline

```bash
nextflow pull rki-mf1/CoVpipe2
```

### Use a certain release

We recommend to use a stable release of the pipeline:

```bash
nextflow pull rki-mf1/CoVpipe2 -r <RELEASE>
```

## Quick run examples

### Example 1:
```bash
nextflow run rki-mf1/CoVpipe2 \
      --reference 'sars-cov-2' \
      --fastq my_samples.csv --list \
      --kraken \
      --cores 4 --max_cores 8
```
- Read input from sample sheet
- Perform taxonomic classification to remove not SARS-CoV-2 reads
- Local execution with maximal 8 cores in total and conda

### Example 2:
```bash
nextflow run rki-mf1/CoVpipe2 \
      --reference 'sars-cov-2' \
      --fastq '*R{1,2}.fastq.gz' \
      --adapter /path/to/repo/data/adapters/NexteraTransposase.fasta \
      --primer_version V4.1 \
      -profile slurm,singularity
```

- Remove adapters
- Clip primer (ARTIC version V4.1)
- Execution on a SLURM system with Singularity

### Example sample sheet

`CoVpipe2` accepts a sample sheet in `CSV` format as input and should look like this:

```
sample,fastq_1,fastq_2
sample1,/path/to/reads/id1_1.fastq.gz,/path/to/reads/id1_2.fastq.gz
sample2,/path/to/reads/id2_1.fastq.gz,/path/to/reads/id2_2.fastq.gz
sample3,/path/to/reads/id3_1.fastq.gz,/path/to/reads/id3_2.fastq.gz
sample4,/path/to/reads/id4_1.fastq.gz,/path/to/reads/id4_2.fastq.gz
```

The header is required. Pay attention the set unique sample names!

## Manual / help

<details><summary>click here to see the complete help message</summary>

```
Robert Koch Institute, MF1 Bioinformatics

    Workflow: CoVpipe2

    Usage examples:
    nextflow run CoVpipe2.nf --fastq '*R{1,2}.fastq.gz' --reference 'sars-cov-2' --cores 4 --max_cores 8
    or
    nextflow run rki-mf1/CoVpipe2 -r <version> --fastq '*R{1,2}.fastq.gz' --ref_genome ref.fasta --cores 4 --max_cores 8

    Reference, required:
    --reference              Currently supported: 'sars-cov-2' (MN908947.3)
    OR
    --ref_genome             Reference FASTA file.
    --ref_annotation         Reference GFF file.

    Illumina read data, required:
    --fastq                  One fastq file 'sample{1,2}.fastq' or multiple: '*.fastq.gz' or '*/*.fastq.gz'
                                 If --list, a csv file
                                 If --dir, a directory containing fastq files

    Optional input settings:
    --list                   Activates csv input for --fastq [default: false]
                                 style and header of the csv is: sample,fastq_1,fastq_2
    OR
    --dir                    Read fastq files from directory [default: false]
    --mode                   Switch between 'paired'- and 'single'-end FASTQ; 'single' is experimental [default: paired]
    --run_id                 Run ID [default: ]


    Adapter clipping:
     --adapter               Define the path of a FASTA file containing the adapter sequences to be clipped. [default: false]

    Trimming and QC:
    --fastp_additional_parameters      Additional parameters for FeatureCounts [default: --qualified_quality_phred 20 --length_required 50]
    
    Taxonomic read filter:
    --kraken                 Activate taxonomic read filtering to exclude reads not classified with specific taxonomic ID (see --taxid) [default: false]
                                 A pre-processed kraken2 database will be automatically downloaded from 
                                 https://zenodo.org/record/3854856 and stored locally.
    --taxid                  Taxonomic ID used together with the kraken2 database for read filtering [default: 2697049]

    Primer detection: 
    --primer_bedpe           Provide the path to the primer BEDPE file. [default: false]
                                 TAB-delimited text file containing at least 6 fields, see here:
                                 https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format
    OR
    --primer_bed             Provide the path to the primer BED file. [default: false]
    OR
    --primer_version         Provide a primer version. Currently supported ARTIC versions: V1, V2, V3, V4, V4.1 [default: false]

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
                                 no strand bias. Set to -1 to disable the filter. [default: 60]
    --var_qual               Minimal variant call quality. Freebayes produces a general judgement of the 
                                 variant call. [default: 10]

    Consensus generation:
    --cns_min_cov            Minimum number of reads required so that the respective position in the consensus sequence 
                                 is NOT hard masked. [default: 20]
    --cns_gt_adjust          Minimum fraction of reads supporting a variant which leads to an explicit call of this 
                                 variant (genotype adjustment). The value has to be greater than 0.5 but not greater than 1. 
                                 To turn genotype adjustment off, set the value to 0. [default: 0.9]

    Updated for linage assignment and mutation calling:
    --update                   Update pangolin and nextclade [default: false]
                                  Depending on the chosen profile either the conda environment (profiles 'standard', 'conda', 'mamba') 
                                  or the container (profiles 'docker', 'singularity') is updated.
    --pangolin_docker_default  Default container tag for pangolin [default: rkimf1/pangolin:3.1.20--3bb06db]
    --nextclade_docker_default Default container tag for nextclade [default: rkimf1/nextclade:1.10.2--1764691]

    Computing options:
    --cores                  Max cores per process for local use [default: 4]
    --max_cores              Max cores used on the machine for local use [default: 12]
    --memory                 Max memory in GB for local use [default: 12]

    Output options:
    --output                 Name of the result folder [default: results]
    --publish_dir_mode       Mode of output publishing: 'copy', 'symlink' [default: copy]
                                 With 'symlink' results are lost when removing the work directory.

    Caching:
    --databases              Location for auto-download data like databases [default: nextflow-autodownload-databases]
    --conda_cache_dir        Location for storing the conda environments [default: conda]
    --singularity_cache_dir  Location for storing the singularity images [default: singularity]
    
    Execution/Engine profiles:
    The pipeline supports profiles to run via different Executers and Engines e.g.: -profile local,conda
    
    Executer (choose one):
      local
      slurm
    
    Engines (choose one):
      conda
      mamba
      docker
      singularity

    Misc:
      cluster                Loads resource configs more suitable for cluster execution.
                             Has to be combine with an engine and an executor.
    

    Per default: -profile local,conda is executed.
```

</details>

## Workflow

Workflow overview:
![workflow](/data/figures/covpipe2_steps.png)
<sub><sub>Components originally designed by James A. Fellows Yates & nf-core under a CC0 license (public domain)</sub></sub>

<details><summary>More detailed overview with process names:</summary>

![workflow](/data/figures/covpipe2_processes.png)
<sub><sub>Components originally designed by James A. Fellows Yates & nf-core under a CC0 license (public domain)</sub></sub>

</details>

<details><summary>Even more detailed overview with process names and parameters:</summary>

![workflow](/data/figures/covpipe2_processes_params.png)
<sub><sub>Components originally designed by James A. Fellows Yates & nf-core under a CC0 license (public domain)</sub></sub>

</details>

## Acknowledgement, props and inspiration

- [ncov_minipipe aka CoVpipe](https://gitlab.com/RKIBioinformaticsPipelines/ncov_minipipe)
- [poreCov](https://github.com/replikation/poreCov)
- [nf-core](https://nf-co.re/pipelines)