/* Comments:
This is the "autodownload" process for the kraken database.
It features auto storage via cloud (line 2) or local storage (line 3)

We use "storeDir" in local mode - because its like a "persistant" process output,
this is automatically always checked by nextflow. Unfortunately, it's not usable for the cloud.
We use the "checkifexists" instead. (see main.nf code)
*/

process kraken_db {
    // label just some environment or container where we are sure that wget and tar are available
    label 'dos2unix'

    if (params.cloudProcess) { publishDir "${params.databases}/kraken2", mode: 'copy', pattern: "GRCh38.p13_SC2_2022-03-01" }
    else { storeDir "${params.databases}/kraken2" }  

    output:
    path("GRCh38.p13_SC2_2022-03-01", type: 'dir')

    script:
    """
    wget https://zenodo.org/record/6333909/files/GRCh38.p13_SC2_2022-03-01.tar.gz?download=1 -O GRCh38.p13_SC2_2022-03-01.tar.gz
    tar zxvf GRCh38.p13_SC2_2022-03-01.tar.gz
    """
    stub:
    """
    touch GRCh38.p13_SC2_2022-03-01
    """
}

process kraken {
    label 'kraken'

    publishDir "${params.output}/${params.read_dir}/${name}/kraken", mode: params.publish_dir_mode, pattern: "*.txt"

    input:
    tuple val(name), path(reads)
    path(db)

    output:
    tuple val(name), file("${name}.classified.R*.fastq"),     emit: fastq
    tuple val(name), file("${name}.kraken.out.txt"),                  emit: kraken_output
    tuple val(name), file("${name}.kraken.report.txt"),               emit: kraken_report

    script:
    set_paired = params.mode == 'paired' ? '--paired' : ''
    set_out_name = params.mode == 'paired' ? '#' : ''
    """
    kraken2 \
        --threads $task.cpus \
        --db ${db} \
        ${set_paired} \
        --classified-out ${name}.classified.R${set_out_name}.fastq \
        --output ${name}.kraken.out.txt \
        --report ${name}.kraken.report.txt \
        --gzip-compressed \
        ${reads}
    """
    stub:
    """
    touch ${name}.classified.R_{1,2}.fastq ${name}.classified.R.fastq ${name}.kraken.out.txt ${name}.kraken.report.txt
    """
}

process filter_virus_reads {
    // label just some environment or container where we are sure that wget and tar are available
    label 'dos2unix'

    publishDir "${params.output}/${params.read_dir}/${name}/kraken/classified", mode: params.publish_dir_mode

    input:
    tuple val(name), path(reads)

    output:
    tuple val(name), file("${name}.classified.R*.fastq.gz"),    emit: fastq
    tuple val(name), file("${name}.extract.log"),                   emit: kraken_report

    script:
    if ( params.mode == 'paired' )
        """
        set +o pipefail
        zgrep -A3 'kraken:taxid|${params.taxid}' ${reads[0]} | sed -e 's/^--\$//' | sed '/^\\s*\$/d' | gzip 1> ${name}.classified.R1.fastq.gz 2>> ${name}.extract.log
        zgrep -A3 'kraken:taxid|${params.taxid}' ${reads[1]} | sed -e 's/^--\$//' | sed '/^\\s*\$/d' | gzip 1> ${name}.classified.R2.fastq.gz 2>> ${name}.extract.log
        """
    else if ( params.mode == 'single' )
        """
        set +o pipefail
        zgrep -A3 'kraken:taxid|${params.taxid}' ${reads} | sed -e 's/^--\$//' | sed '/^\\s*\$/d' | gzip 1> ${name}.classified.R.fastq.gz 2>> ${name}.extract.log
        """
    stub:
    """
    touch ${name}.classified.R{1,2}.fastq.gz ${name}.extract.log
    """
}
