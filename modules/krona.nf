process krona {
    label 'krona'   
    publishDir "${params.output}/${params.read_dir}/${name}/kraken", mode: params.publish_dir_mode
    errorStrategy 'ignore'

    input:
        tuple val(name), path(kreport)
  	output:
    	tuple val(name), file("${name}_krona.html")
    
  	script:
    """
    cat ${kreport} | cut -f 3,5 > file.krona
    # ktUpdateTaxonomy.sh
    ktImportTaxonomy file.krona -m 1
    mv *.html ${name}_krona.html
    """
    stub:
    """
    touch ${name}_krona.html
    """
}