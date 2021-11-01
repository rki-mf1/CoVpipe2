process compress_reads {
  label 'pigz'

  input:
  tuple val(name), path(reads)
  
  output:
  tuple val(name), path("*.fastq.gz")

  script:
  """
  pigz -f -p ${task.cpus} ${reads}
  """
}

process bgzip_compress {
  label 'samtools'

  publishDir "${params.publish_dir}/${name}", mode: params.publish_dir_mode

  input:
  tuple val(name), path(file)

  output:
  tuple val(name), path("${file}.gz")

  script:
  """
  bgzip -@ ${task.cpus} ${file}
  """
}

process adapt_consensus_header {

  input:
  tuple val(name), path(fasta)

  output:
  tuple val(name), path("${name}.iupac_consensus.fasta")

  script:
  """
    VERSION='unknown_version'
    if [ ${workflow.revision} != 'null' ]; then
      VERSION=${workflow.revision}
    fi
    sed "1 s/.*/>${name}_iupac_consensus_\${VERSION}/" ${fasta} 1> ${name}.iupac_consensus.fasta
  """
  }


process mask_iupac {
  input:
  tuple val(name), path(fasta)

  output:
  tuple val(name), path("${name}.masked_consensus.fasta")

  script:
    """
      VERSION='unknown_version'
      if [ ${workflow.revision} != 'null' ]; then
        VERSION=${workflow.revision}
      fi
      echo ">${name}_masked_consensus_\${VERSION}\" 1> ${name}.masked_consensus.fasta 2> ${name}_createMaskedConsensus.log
      tail -n +2 ${fasta} | tr "RYSWKMBDHVN" "N" 1>> ${name}.masked_consensus.fasta 2>> ${name}_createMaskedConsensus.log
    """
}