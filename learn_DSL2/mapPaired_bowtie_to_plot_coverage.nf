#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/* Usage
     nextflow run word_count.nf --reference <fasta file> --inputPair <paired fastq files>
 */

params.reference = "${baseDir}/bowtieConsensTestFiles/eng_live_atten_poliovirus/MZ245455.1"
reference_name = file(params.reference).name
reference_idx = "${reference_name}.fai"
reference_path = file(params.reference).parent


params.inputPair = "${baseDir}/bowtieConsensTestFiles/eng_live_atten_poliovirus/polio_sample_3_screened_trim_R?_001.fastq.gz"
params.output = "${baseDir}/../output/"

process LOOKSY  // for debugging and sanity checking
  {

    input:
    tuple val(sample_name), path(reads)
    path(reference)

    output:
    stdout

    script:
    """
     printf '${reference} first line is: -> '
     head -1 ${reference}
     if [[ "${reads[0]}" =~ .*gz\$ ]] && [[ "${reads[1]}" =~ .*gz\$ ]];
     then
        printf '${reads[0]} first line is: ->  '
        zcat ${reads[0]} | head -1
        printf '${reads[1]} first line is: ->  '
        zcat ${reads[1]} | head -1
     else
        printf '${reads[0]} first line is ->  '
        head -1 ${reads[0]}
        printf '${reads[1]} first line is: ->  '
        head -1 ${reads[1]} 
     fi
    """
  
  }


process bowtie2map {

    publishDir "${baseDir}/../intermediate/", mode: 'copy'
    
    input:
    tuple val(sample_name), path(reads)
    path reference 

    output:
    tuple val(sample_name), path("${sample_name}.sam")    

    script:
    """
    bowtie2 --no-unal --no-mixed -x "${reference}"/"${reference_name}" -1 "${reads[0]}" -2 "${reads[1]}" > "${sample_name}.sam"
    """
}


process bowtie2map_singularity {

    publishDir "${baseDir}/../intermediate/", mode: 'copy'
    
    input:
    tuple val(sample_name), path(reads)
    path reference 

    output:
    tuple val(sample_name), path("${sample_name}.sam")    

    script:
    """
    singularity exec "${baseDir}/../"my_bowtie2.sif bowtie2 --no-unal --no-mixed -x "${reference}"/"${reference_name}" -1 "${reads[0]}" -2 "${reads[1]}" > "${sample_name}.sam"
    """
}

process sam2bam {

    publishDir "${baseDir}/../intermediate/", mode: 'copy'
    
    input:
    tuple val(sample_name), path("${sample_name}.sam")
    
    output:
    tuple val(sample_name), path("${sample_name}.bam")

    script:
    """
    samtools view "${sample_name}".sam -o "${sample_name}".bam
    """

}

process sortBam {

    publishDir "${baseDir}/../output/", mode: 'copy'

    input:
    tuple val(sample_name), path("${sample_name}.bam")
    
    output:
    tuple val(sample_name), path("${sample_name}.sorted.bam")
    
    script:
    """
    samtools sort "${sample_name}".bam -o "${sample_name}".sorted.bam
    """
}


process indexBam {
    
    publishDir "${baseDir}/../output/", mode: 'copy'
    
    input:
    tuple val(sample_name), path("${sample_name}.sorted.bam")
    
    output:
    tuple val(sample_name), path("${sample_name}.sorted.bam.bai")
    
    script:
    """
    samtools index -b "${sample_name}".sorted.bam "${sample_name}".sorted.bam.bai
    """

}

process makeGenomeCov {

   publishDir "${baseDir}/../output/", mode: 'copy'
   
   input:
   tuple val(sample_name), path("${sample_name}.sorted.bam")
   path reference
   
   output:
   tuple val(sample_name), path("${sample_name}.bedGraph")
   
   script:
   """
   genomeCoverageBed -ibam "${sample_name}".sorted.bam -d -g "${reference}"/"${reference_name}".sizes > "${sample_name}".bedGraph
   """
   
 }


process callPerl {
  
  publishDir "${baseDir}/../output/", mode: 'copy'

  input:
  tuple val(sample_name), path("${sample_name}.bedGraph")

  output:
  stdout

  script:
  """
  perl ${baseDir}/../Perl_scripts/coverageStatsFromBedGraph.pl "${sample_name}".bedGraph
  """
  
}


workflow {
    
    Channel
        .fromFilePairs(params.inputPair, checkIfExists: true)
        .set { read_pairs_ch }


   // LOOKSY(read_pairs_ch, params.reference) | view
    
   //  mapResults = bowtie2map(read_pairs_ch, reference_path)

   mapResults = bowtie2map(read_pairs_ch, reference_path) 

   mapResults.view { "Bowtie2 Results: ${it}" }

   bamResults = sam2bam(mapResults)

   sortedBam = sortBam(bamResults)

   sortedBam.view { "sortBam results: ${it}" }
   
   indexedBam = indexBam(sortedBam)

   bedGr = makeGenomeCov(sortedBam, reference_path)

   bedGr.view { "Coverage Plot: ${it}" }

   callPerl(bedGr).view()

}
