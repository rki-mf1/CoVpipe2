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