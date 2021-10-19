process index_bwa {
    label 'bwa'  
    
    publishDir "${params.output}", mode: 'copy', pattern: "${name}_bwa_index.tar.gz"
    
    input:
    tuple val(name), file(reference)
    
    output:
    tuple val(name), file("${name}_bwa_index.tar.gz"), emit: index

    script:
    """
    bwa index ${reference} &> /dev/null
    tar -zcvf ${name}_bwa_index.tar.gz ${reference}.*
    """
}

