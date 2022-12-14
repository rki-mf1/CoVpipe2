process fastp {
    label 'fastp'

    publishDir "${params.output}/${params.read_dir}/${name}/trimming", mode: params.publish_dir_mode

    input:
        tuple val(name), path(reads)
        path(adapter)

    output:
        tuple val(name), path("${name}.fastp.R{1,2}.fastq.gz"),       emit: reads
        tuple val(name), path("${name}.fastp.json"),            emit: json
        tuple val(name), path("${name}.fastp.html"),            emit: html
    
    script:
    set_adapters = adapter.getName() == 'NO_ADAPTERS' ? '' : "--adapter_fasta ${adapter}"
    set_paired_reads = params.mode == 'single' ? '' : "--in2 ${reads[1]} --out2 ${name}.fastp.R2.fastq.gz --unpaired1 ${name}.SE.R1.fastq.gz --unpaired2 ${name}.SE.R2.fastq.gz"
    """
    fastp \
        --in1 ${reads[0]} \
        --out1 ${name}.fastp.R1.fastq.gz \
        ${set_paired_reads} \
        --json ${name}.fastp.json \
        --html ${name}.fastp.html \
        ${set_adapters} \
        --low_complexity_filter \
        --overrepresentation_analysis \
        --thread $task.cpus \
        ${params.fastp_additional_parameters}
    """
    stub:
    """
    touch ${name}.fastp.R{1,2}.fastq.gz ${name}.fastp.json ${name}.fastp.html
    """
}
