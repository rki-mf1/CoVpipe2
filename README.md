# CoVpipe2

CoVpipe2 is a Nextflow pipeline for reference-based genome reconstruction of SARS-CoV-2 from NGS data. In principle it can be used also for other viruses.

<details><summary>Table of contents</summary>

- [CoVpipe2](#covpipe2)
  - [Quick installation](#quick-installation)
    - [Call help](#call-help)
    - [Test run](#test-run)
    - [Update the pipeline](#update-the-pipeline)
    - [Use a certain release](#use-a-certain-release)
  - [Quick run examples](#quick-run-examples)
    - [Example 1:](#example-1)
    - [Example 2:](#example-2)
    - [Example sample sheet](#example-sample-sheet)
  - [Manual](#manual)
  - [Changes to CoVpipe](#changes-to-covpipe)
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

### Test run

Validate your installation with a test run:

```bash
# for a Conda installation
nextflow run rki-mf1/CoVpipe2 -profile local,conda,test --cores 4 --max_cores 8

# for a Singularity installation
nextflow run rki-mf1/CoVpipe2 -profile local,singularity,test --cores 4 --max_cores 8

# for a Docker installation
nextflow run rki-mf1/CoVpipe2 -profile local,docker,test --cores 4 --max_cores 8
```

For more configuration options, see [here](#manual).

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

## Manual

<details><summary>click here to see the complete help message</summary>

```
Robert Koch Institute, MF1 Bioinformatics

    Workflow: CoVpipe2

    Usage examples:
    nextflow run CoVpipe2.nf --fastq '*R{1,2}.fastq.gz' --cores 4 --max_cores 8
    or
    nextflow run rki-mf1/CoVpipe2 -r <version> --fastq '*R{1,2}.fastq.gz' --ref_genome ref.fasta --cores 4 --max_cores 8

    Reference, required:
    --reference              Currently supported: 'sars-cov-2' (MN908947.3) [default: sars-cov-2]
    OR
    --ref_genome             Reference FASTA file.
    --ref_annotation         Reference GFF file.

    Illumina read data, required:
    --fastq                  e.g.: 'sample{1,2}.fastq' or '*.fastq.gz' or '*/*.fastq.gz'

    Optional input settings:
    --list                   This flag activates csv input for --fastq [default: false]
                                 style and header of the csv is: sample,fastq_1,fastq_2
    --mode                   Switch between 'paired'- and 'single'-end FASTQ; 'single' is experimental [default: paired]
    --run_id                 Run ID [default: ]

    Adapter clipping:
     --adapter               Define the path of a FASTA file containing the adapter sequences to be clipped. [default: false]

    Trimming and QC:
    --fastp_additional_parameters      Additional parameters for fastp [default: --qualified_quality_phred 20 --length_required 50]
                                           For shorter/longer amplicon length than 156 nt, adjust --length_required
    
    Taxonomic read filter:
    --kraken                 Activate taxonomic read filtering to exclude reads not classified with specific taxonomic ID (see --taxid) [default: false]
                                 A pre-processed kraken2 database will be automatically downloaded from 
                                 https://zenodo.org/record/3854856 and stored locally.
    --kraken_db_custom       Path to a custom Kraken2 database. [default: ]
    --taxid                  Taxonomic ID used together with the kraken2 database for read filtering [default: 2697049]

    Linage detection on read level with LCS:
    Uses this fork https://github.com/rki-mf1/LCS of https://github.com/rvalieris/LCS
    --read_linage            Linage detection on read level [default: false]
    --lcs_ucsc_version       Create marker table based on a specific UCSC SARS-CoV-2 tree (e.g. '2022-05-01'). Use 'predefined' 
                                 to use the marker table from the repo (most probably not up-to-date) [default: predefined]
                                 See https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2 for available trees.
    --lcs_ucsc_predefined    If '--lcs_ucsc_version 'predefined'', select pre-calculated UCSC table [default: 2022-01-31]
                                 See https://github.com/rki-mf1/LCS/tree/master/data/pre-generated-marker-tables
    --lcs_ucsc_update        Use latest UCSC SARS-CoV-2 tree for marker table update. Overwrites --lcs_ucsc_version [default: false]
                                 Automatically checks https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/public-latest.version.txt
    --lcs_ucsc_downsampling  Downsample sequences when updating marker table to save resources. Use 'None' to turn off [default: 10000]
                                 Attention! Updating without downsampling needs a lot of resources in terms of memory and might fail.
                                 Consider downsampling or increase the memory for this process.
    --lcs_variant_groups     Provide path to custom variant groups table (TSV) for marker table update. Use 'default' for predefined groups from repo
                                 (https://github.com/rki-mf1/LCS/blob/master/data/variant_groups.tsv) [default: default]
    --lcs_cutoff             Plot linages above this threshold [default: 0.03]

    Mapping: 
    --isize_filter           Insert size threshold for mapping. All BAM file entries with an insert size above this threshold 
                                 are filtered out. Deactivated by default. [default: false]

    Primer detection: 
    --bamclipper_additional_parameters      Additional parameters for BAMClipper [default: false]
                                                Use -u INT and -d INT to adjust the primer detection window of BAMClipper: extend upstream (-u) or 
                                                downstream (-d) from the 5' most nt of primer [default from BAMClipper: -u 1 -d 5]
    --primer_bedpe           Provide the path to the primer BEDPE file. [default: false]
                                 TAB-delimited text file containing at least 6 fields, see here:
                                 https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format
    OR
    --primer_bed             Provide the path to the primer BED file. A BEDPE file will be generated automatically.
                                 The name of each entry has to match this pattern: primerID[_LEFT|_RIGHT]_ampliconID [default: false]
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
    --var_sap                Maximal strand balance probability for the alternate allele (SAP). The SAP is the Phred-scaled 
                                 probability that there is strand bias at the respective site. A value near 0 indicates little or 
                                 no strand bias. Amplicon data usually has a high, WGS data a low bias. [default: false]
                                 Disable (default) for amplicon sequencing; for WGS GATK recommends 60
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
    The pipeline supports profiles to run via different Executors and Engines e.g.: -profile local,conda
    
    Executor (choose one):
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

    Test profile:
    Test the pipeline with a small test dataset:
    nextflow run rki-mf1/CoVpipe2 -profile executor,engine,test
```

</details>

## Changes to [CoVpipe](https://gitlab.com/RKIBioinformaticsPipelines/ncov_minipipe)

- Workflow management framework: `snakemake` -> `Nextflow`
- Docker/Singularity and conda support for each step
- Container/conda updated feature for `pangolin` and `nextclade`
- HPC/slurm profile provided
- Fixes:
  - Subtract only deletions from low coverage mask for consensus generation
- New features:
  - `nexclade` (mutation calling, clade assignment)
  - `LCS` (linage decomposition)
  - Restructured report
  - `krona` plots (visualization of `Kraken2` output)
  - `president` (genome quality control)
- Version update (status CoVpipe2 v0.2.1):
  - `bcftools`: 1.11 -> 1.14
    - Note: https://github.com/samtools/bcftools/issues/1708
  - `liftoff`: 1.5.2 -> 1.6.2
  - `kraken2`: 2.1.0 -> 2.1.2
  - `freebayes`: 1.3.2 -> 1.3.6
  - `fastp`: 0.20.1 -> 0.23.2
  - `bedtools`: 2.29.2 -> 2.30.0

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