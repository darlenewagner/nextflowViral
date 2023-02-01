#! /apps/x86_64/nextflow/21.04.3 nextflow

// This pipeline takes a Bowtie2-indexed reference and an arbitrary number of paired-end
// .fastq files to create reference-based assembly .fasta files in params.consensdir
// Sorted bam pileups and VCF files are also created and stored in params.intermdir
// For reference-based assemblies, bases with less than 1X or 25X coverage are masked with 'N'
// Indels and complex polymorphisms are excluded from the VCF files
// Prerequisites: nextflow, bowtie2, samtools, freebayes, and genomeCoverageBed

Channel
  .fromFilePairs(params.querydir)
  .set { samples_ch }
  
process mapSAM {

  publishDir "$params.intermdir", mode: 'copy'  
  
  input:
  tuple val(sampleId), path(reads) from samples_ch
  
  output:
  tuple sampleId, path("${sampleId}.sam") into samOut
  
  script:
  def (read1, read2) = reads
  if( params.strictPairing == true )
    """
    /apps/x86_64/bowtie2/bowtie2-2.3.5.1/bowtie2-2.3.5.1-linux-x86_64/bowtie2 --no-unal --no-mixed -x $params.reference -1 "${read1}" -2 "${read2}" -S "${sampleId}.sam"
    """
  
  else
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
  tuple sampleId, path("${sampleId}.sorted.bam") into sortedOut0

  script:
  """
  samtools sort "${sampleId}".bam -o "${sampleId}".sorted.bam  
  """
}

process makeGenomeCov {
  
  input:
  tuple sampleId, path("${sampleId}.sorted.bam") from sortedOut0

  output:
  tuple sampleId, path("${sampleId}.bedGraph") into bedOut
  tuple sampleId, path("${sampleId}.bedGraph") into bedOut2

  script:
  """
  genomeCoverageBed -ibam "${sampleId}".sorted.bam -d -g "${params.reference}".sizes > "${sampleId}".bedGraph
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

process makeVCF {

  publishDir "$params.intermdir", mode: 'copy'
    
  input:
  tuple sampleId, path("${sampleId}.sorted.bam") from sortedOut2
  
  output:
  tuple sampleId, path("${sampleId}.vcf") into bVCF
//  tuple sampleId, path("${sampleId}.mpileup") into bVCF2
  
  script:
  if( params.freebayes == true )
  """
  freebayes --min-coverage 5 --no-indels --no-mnps --no-complex -f "${params.reference}".fasta "${sampleId}".sorted.bam > "${sampleId}".vcf
  """
  
  else
  """
  bcftools mpileup -d 35000 -Ob -f "${params.reference}".fasta -Q 20 -q 20 --annotate FORMAT/AD,FORMAT/DP "${sampleId}".sorted.bam | bcftools call -mv --ploidy 1 -Ov --output "${sampleId}".vcf 
  """
}

process prepVCF {

  publishDir "$params.intermdir", mode: 'copy'
  
  input:
  tuple sampleId, path("${sampleId}.vcf") from bVCF
  
  output:
  tuple sampleId, path("${sampleId}.vcf.gz") into bgzVCF
  tuple sampleId, path("${sampleId}.vcf.gz") into bgzVCF2
  tuple sampleId, path("${sampleId}.vcf.gz.csi") into idxVCF
  
  script:
  """
  bgzip -c "${sampleId}".vcf > "${sampleId}".vcf.gz  
  bcftools index "${sampleId}".vcf.gz -o "${sampleId}".vcf.gz.csi
  """

}

process makeBcfConsensus {
  
  input:
  tuple sampleId, path("${sampleId}.vcf.gz") from bgzVCF
  tuple sampleId, path("${sampleId}.vcf.gz.csi") from idxVCF
  
  output:
  tuple sampleId, path("${sampleId}.fasta") into FASTA
  tuple sampleId, path("${sampleId}.fasta") into FASTA2
  
  script:
  """
  cat "${params.reference}".fasta | bcftools consensus "${sampleId}".vcf.gz > "${sampleId}".fasta
  """
}

process makeCoverageMask {

  input:
  tuple sampleId, path("${sampleId}.bedGraph") from bedOut

  output:
  tuple sampleId, path("${sampleId}.bedGraph.Zero.nnn") into zeroNNN
  
  script:
  """
  awk '{ if (\$3 < 1) {print \$1"\t" \$2"\t" \$3"\t" "N" }}' "${sampleId}".bedGraph > "${sampleId}".bedGraph.Zero.nnn 
  """
}

process makeCoverageMask2 {
  
  input:
  tuple sampleId, path("${sampleId}.bedGraph") from bedOut2
  
  output:
  tuple sampleId, path("${sampleId}.bedGraph.Cov25X.nnn") into twentyFiveNNN
  
  script:
  """
  awk '{ if (\$3 < 25) {print \$1"\t" \$2"\t" \$3"\t" "N" }}' "${sampleId}".bedGraph > "${sampleId}".bedGraph.Cov25X.nnn 
  """
}

process maskWithNs {

  input:
  tuple sampleId, path("${sampleId}.bedGraph.Zero.nnn") from zeroNNN
  tuple sampleId, path("${sampleId}.fasta") from FASTA

  output:
  tuple sampleId, path("${sampleId}.atLeast1X.fa") into atLeast1X
  
  script:
  """
  seqtk mutfa "${sampleId}".fasta "${sampleId}".bedGraph.Zero.nnn > "${sampleId}".atLeast1X.fa   
  """
}

process maskWithNs2 {
  
  input:
  tuple sampleId, path("${sampleId}.bedGraph.Cov25X.nnn") from twentyFiveNNN  
  tuple sampleId, path("${sampleId}.fasta") from FASTA2

  output:
  tuple sampleId, path("${sampleId}.atLeast25X.fa") into atLeast25X  
  
  script:
  """
  seqtk mutfa "${sampleId}".fasta "${sampleId}".bedGraph.Cov25X.nnn > "${sampleId}".atLeast25X.fa     
  """
}


process renameHeader {

  publishDir "$params.consensdir", mode: 'copy'
  
  input:
  tuple sampleId, path("${sampleId}.atLeast1X.fa") from atLeast1X

  output:
  tuple sampleId, path("${sampleId}.atLeast1X.fasta") into atLeast1XX

  script:
  """
  sed -i "1s/.*/>"${sampleId}"/" "${sampleId}".atLeast1X.fa
  cat "${sampleId}".atLeast1X.fa > "${sampleId}".atLeast1X.fasta
  """
}


process renameHeader2 {

  publishDir "$params.consensdir", mode: 'copy'
  
  input:
  tuple sampleId, path("${sampleId}.atLeast25X.fa") from atLeast25X  
  
  output:
  tuple sampleId, path("${sampleId}.atLeast25X.fasta") into atLeast25XX  
  
  script:
  """
  sed -i "1s/.*/>"${sampleId}"/" "${sampleId}".atLeast25X.fa
  cat "${sampleId}".atLeast25X.fa > "${sampleId}".atLeast25X.fasta 
  """
}


process queryVCF {
  
  publishDir "$params.consensdir", mode: 'copy'
  
  input:
  tuple sampleId, path("${sampleId}.vcf.gz") from bgzVCF2
  
  output:
  tuple sampleId, path("${sampleId}.snp.tsv") into tabSNP
  
  script:
  """
  bcftools query -f '''%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t[%AD]\t[%DP]\n''' "${sampleId}".vcf.gz | perl $PWD/bcftoolsQuery.pl > "${sampleId}".snp.tsv
  """
}
