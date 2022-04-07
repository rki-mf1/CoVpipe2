process krona {
    label 'krona'
    publishDir "${params.output}/${params.read_dir}/${name}/kraken", mode: params.publish_dir_mode
    if ( workflow.profile.contains('mamba') ||  workflow.profile.contains('conda') ){ errorStrategy 'ignore' }

    input:
        tuple val(name), path(kreport)
        val(tax_status)
  	output:
    	tuple val(name), file("${name}_krona.html")
    
  	script:
    """
    cat ${kreport} | cut -f 3,5 > file.krona
    ktImportTaxonomy file.krona -m 1
    mv *.html ${name}_krona.html
    """
    stub:
    """
    touch ${name}_krona.html
    """
}

process krona_taxonomy_update {
    label 'krona'
    executor 'local'

    output:
        val('updated')

    script:
    """
    ktUpdateTaxonomy.sh
    """
}