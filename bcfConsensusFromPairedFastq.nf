#! /apps/x86_64/nextflow/21.04.3 nextflow

// This pipeline takes a SMALT-indexed reference and an arbitrary number of paired-end
// .fastq files to create a vcf file of SNPs. Indels and complex polymorphisms are excluded 
// from the .vcf requirements: nextflow, smalt, samtools, and freebayes

Channel
  .fromFilePairs('paired_files/*trim_dedup_R{1,2}.fastq')
  .set { samples_ch }
  
process mapSAM {

//  publishDir "/scicomp/home-pure/ydn3/nextflow_for_read_mapping/mapOut/", mode: 'copy'  
  
  input:
  tuple val(sampleId), path(reads) from samples_ch
  
  output:
  tuple sampleId, path("${sampleId}.sam") into samOut
  
  script:
  def (read1, read2) = reads
  """
  smalt map -o "${sampleId}".sam $params.reference "${read1}" "${read2}" 
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
  
  input:
  tuple sampleId, path("${sampleId}.bam") from bamOut
  
  output:
  tuple sampleId, path("${sampleId}.sorted.bam") into sortedOut
  tuple sampleId, path("${sampleId}.sorted.bam") into sortedOut2

  script:
  """
  samtools sort "${sampleId}".bam -o "${sampleId}".sorted.bam  
  """
}

process indexBam {
  
  input:
  tuple sampleId, path("${sampleId}.sorted.bam") from sortedOut
  
  output:
  tuple sampleId, path("${sampleId}.sorted.bai") into indexOut
  
  script:
  """
  samtools index "${sampleId}".sorted.bam "${sampleId}".sorted.bai
  """
}

process makeVCF {

  publishDir "/scicomp/home-pure/ydn3/nextflow_for_read_mapping/mapOut/", mode: 'copy'
    
  input:
  tuple sampleId, path("${sampleId}.sorted.bam") from sortedOut2
  
  output:
  tuple sampleId, path("${sampleId}.vcf") into bVCF
  
  script:
  """
  freebayes --min-coverage 5 --no-indels --no-mnps --no-complex -f $params.rawReference "${sampleId}".sorted.bam >"${sampleId}".vcf  
  """
}

process prepVCF {
  publishDir "/scicomp/home-pure/ydn3/nextflow_for_read_mapping/mapOut/", mode: 'copy'
  
  input:
  tuple sampleId, path("${sampleId}.vcf") from bVCF
  
  output:
  tuple sampleId, path("${sampleId}.vcf.gz") into bgzVCF
  tuple sampleId, path("${sampleId}.vcf.gz.csi") into idxVCF
  
  script:
  """
  bgzip -c "${sampleId}".vcf > "${sampleId}".vcf.gz  
  bcftools index "${sampleId}".vcf.gz -o "${sampleId}".vcf.gz.csi
  """
}

process makeBcfConsensus {
  publishDir "/scicomp/home-pure/ydn3/nextflow_for_read_mapping/multi_consensus/", mode: 'copy'
  
  input:
  tuple sampleId, path("${sampleId}.vcf.gz") from bgzVCF
  tuple sampleId, path("${sampleId}.vcf.gz.csi") from idxVCF
  
  output:
  tuple sampleId, path("${sampleId}.fasta") into FASTA
  
  script:
  """
  cat $params.rawReference | bcftools consensus "${sampleId}".vcf.gz > "${sampleId}".fasta
  """
}
