#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/* Usage for test data
nextflow run mapSingle_bowtie_to_plot_coverage.nf --reference "$PWD/../bowtieConsensTestFiles/adenovirus_B3/OR777202.1" --inputSingle "$PWD/../bowtieConsensTestFiles/adenovirus_B3/Pool-1_S1_adenovirus_B3_001.fastq.gz" --intermediate "$PWD/intermediate/"
 */

params.reference = "${baseDir}/../bowtieConsensTestFiles/eng_live_atten_poliovirus/MZ245455.1"
reference_name = file(params.reference).name
reference_idx = "${reference_name}.fai"
reference_path = file(params.reference).parent

params.inputSingle = "${baseDir}/../expBowtieConsensTestFiles/pseudoPairs/Pool-1_S1_final_pseudoPairs.fastq"
inputSingle_path = file(params.inputSingle).parent
params.output = "${baseDir}/output/"
params.intermediate = "${baseDir}/intermediate/"

process LOOKSY  // for debugging and sanity checking
  {

    input:
    path(sample_name)
    path(reference)
    
    output:
    stdout
    
    script:
    """
     printf '${reference} first line is: -> '
     head -1 ${reference}
     if [[ "${sample_name}" +~ .*fastq\$ ]];
     then
        printf '${sample_name} first line is: ->  '
        zcat ${sample_name} | head -1
        printf '${sample_name} first line is: ->  '
        cat ${sample_name} | head -1
     fi
    """
  
  }


process bowtie2map {

    publishDir "${baseDir}/intermediate/", mode: 'copy'
    
    input:
    path(sample_name)
    path reference 

    output:
    tuple val(sample_name.baseName), path("${sample_name.baseName}.sam")    

    script:
    """
    bowtie2 --no-unal -x "${reference}"/"${reference_name}" -U "${sample_name}" -S "${sample_name.baseName}.sam"
    """
}


process bowtie2map_singularity {

    tag { sample }
    
    publishDir "${baseDir}/../intermediate/", mode: 'copy'
    
    input:
    path(sample_name)
    path reference 
    
    output:
    tuple val(sample_name.baseName), path("${sample_name.baseName}.sam")    
       
    script:
    """
    singularity exec "${baseDir}/../"my_bowtie2.sif bowtie2 --no-unal --no-mixed -x "${reference}"/"${reference_name}" -U "${sample_name}" > "${sample_name.baseName}.sam"
    """
}


process sam2bam {

    publishDir "${baseDir}/intermediate/", mode: 'copy'
    
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

    publishDir "${baseDir}/../intermediate/", mode: 'copy'

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
    
    publishDir "${baseDir}/../intermediate/", mode: 'copy'
    
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
        .fromPath(params.inputSingle, checkIfExists: true)
        .set { fastq_ch }
   
   
   // LOOKSY(fastq_ch, params.reference) | view
    
   //  mapResults = bowtie2map(fastq_ch, reference_path)
    
    mapResults = bowtie2map(fastq_ch, reference_path) 
    
    mapResults.view { "Bowtie2 Results: ${it}" }
    
    bamResults = sam2bam(mapResults)
    
    sortedBam = sortBam(bamResults)
    
    indexedBam = indexBam(sortedBam)
    
    bedGr = makeGenomeCov(sortedBam, reference_path)

    bedGr.view { "Coverage Plot: ${it}" }

    callPerl(bedGr).view()   
}
