process lcs_sc2 {
    label 'lcs_sc2'
    publishDir "${params.output}/${params.read_dir}/${name}/lineage-proportion-by-reads", mode: 'copy'
    
    input:
    tuple val(name), path(reads)
    
    output:
    tuple val(name), path("${name}.lcs.tsv")
    
    script:
    """
    git clone https://github.com/rvalieris/LCS.git
    mkdir -p LCS/outputs/variants_table
    zcat LCS/data/pre-generated-marker-tables/pango-designation-markers-v1.2.124.tsv.gz > LCS/outputs/variants_table/pango-markers-table.tsv
    zcat LCS/data/pre-generated-marker-tables/ucsc-markers-2022-01-31.tsv.gz > LCS/outputs/variants_table/ucsc-markers-table.tsv
    mkdir -p LCS/data/fastq
    if [ ${params.mode} == 'paired' ]; then
        cp ${reads[0]} LCS/data/fastq/${name}_1.fastq.gz
        cp ${reads[1]} LCS/data/fastq/${name}_2.fastq.gz
    else
        cp ${reads} LCS/data/fastq/${name}.fastq.gz
    fi
    echo ${name} > LCS/data/tags_pool_mypool
    cd LCS
    mem=\$(echo ${task.memory} | cut -d' ' -f1)
    snakemake --config markers=ucsc dataset=mypool --cores ${task.cpus} --resources mem_gb=\$mem --set-threads pool_mutect=${task.cpus}
    cd ..
    cp LCS/outputs/decompose/mypool.out ${name}.lcs.tsv

    rm -rf LCS/data/fastq
    """
    stub:
    """
    touch ${name}.lcs.tsv
    """
}

process lcs_plot {
    label "r"
    publishDir "${params.output}/${params.read_dir}/", mode: 'copy'

    input:
    path(tsv)
    val(cutoff)
    
    output:
    path("*.png")
    
    script:
    """
    lcs_bar_plot.R '${tsv}' ${cutoff}
    """
}