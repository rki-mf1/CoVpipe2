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