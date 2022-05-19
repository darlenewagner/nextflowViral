#! /apps/x86_64/nextflow/21.04.3 nextflow

// This pipeline takes a Bowtie2-indexed reference and an arbitrary number of paired-end
// .fastq files to calculate read coverage and insert lengths in params.statsdir
// Prerequisites: nextflow, bowtie2, and samtools

Channel
  .fromFilePairs(params.querydir)
  .set { samples_ch }
  
process mapSAM {

  publishDir "$params.intermdir", mode: 'copy'  
  
  input:
  tuple val(sampleId), path(reads) from samples_ch
  
  output:
  tuple sampleId, path("${sampleId}.sam") into samOut
  tuple sampleId, path("${sampleId}.sam") into samOut2
  
  script:
  def (read1, read2) = reads
  """
  /apps/x86_64/bowtie2/bowtie2-2.3.5.1/bowtie2-2.3.5.1-linux-x86_64/bowtie2 --no-unal --local -x $params.reference -1 "${read1}" -2 "${read2}" -S "${sampleId}.sam"
  """
  
}

process sam2bam {
  
  input:
  tuple sampleId, path("${sampleId}.sam") from samOut
  
  output:
  tuple sampleId, path("${sampleId}.bam") into bamOut
  
  script:
  """
  samtools view "${sampleId}".sam -o "${sampleId}".bam  
  """
}

process sortBam {
  
  publishDir "$params.intermdir", mode: 'copy'  

  input:
  tuple sampleId, path("${sampleId}.bam") from bamOut
  
  output:
  tuple sampleId, path("${sampleId}.sorted.bam") into sortedOut
  tuple sampleId, path("${sampleId}.sorted.bam") into sortedOut2
//  tuple sampleId, path("${sampleId}.sorted.bam") into sortedOut0

  script:
  """
  samtools sort "${sampleId}".bam -o "${sampleId}".sorted.bam  
  """
}



process indexBam {
  
  publishDir "$params.intermdir", mode: 'copy'  

  input:
  tuple sampleId, path("${sampleId}.sorted.bam") from sortedOut
  
  output:
  tuple sampleId, path("${sampleId}.sorted.bai") into indexOut
  
  script:
  """
  samtools index "${sampleId}".sorted.bam "${sampleId}".sorted.bai
  """
}


process getInserts 
{
  
  publishDir "$params.statsdir", mode: 'copy'  
  
  input:
  tuple sampleId, path("${sampleId}.sam") from samOut2
  
  output:
  tuple sampleId, path("${sampleId}.INSERT.txt") into inserted  
  
  script:
  """
  perl $PWD/readForInserts.pl "${sampleId}".sam > "${sampleId}".INSERT.txt
  """
}

process getCoverage
{
  publishDir "$params.statsdir", mode: 'copy'  
  
  input:
  tuple sampleId, path("${sampleId}.sorted.bam") from sortedOut2
  
  output:
  tuple sampleId, path("${sampleId}.COVER.txt") into coverage  
  
  script:
  """
  samtools depth -a "${sampleId}".sorted.bam | awk '{c++;s+=\$3}END{print s/c}' > "${sampleId}".COVER.txt
  """
}