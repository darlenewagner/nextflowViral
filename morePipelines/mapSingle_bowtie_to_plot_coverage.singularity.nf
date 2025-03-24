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


process bowtie2map_singularity {

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/bowtie2:2.5.4--he96a11b_5' :
    'quay.io/biocontainers/bowtie2:2.5.4--he96a11b_5' }"
    
    publishDir "${baseDir}/../intermediate/", mode: 'copy'
    
    input:
    path(sample_name)
    path reference 

    output:
    tuple val(sample_name.baseName), path("${sample_name.baseName}.single.sam")    

    script:
    """
    bowtie2 --no-unal -x "${reference}"/"${reference_name}" -U "${sample_name}" -S "${sample_name.baseName}.single.sam"
    """
}


process sam2bam_singularity {

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/samtools:1.9--h91753b0_8' :
    'quay.io/biocontainers/samtools:1.9--h91753b0_8' }"

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

process sortBam_singularity {

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/samtools:1.9--h91753b0_8' :
    'quay.io/biocontainers/samtools:1.9--h91753b0_8' }"

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


process indexBam_singularity {
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/samtools:1.9--h91753b0_8' :
    'quay.io/biocontainers/samtools:1.9--h91753b0_8' }"

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

process makeGenomeCov_singularity {

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/bedtools:2.27.1--he941832_2' :
    'quay.io/biocontainers/bedtools:2.27.1--he941832_2' }"   

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
   
   
    mapResults = bowtie2map_singularity(fastq_ch, reference_path) 
    
    mapResults.view { "Bowtie2 Results: ${it}" }
    
    bamResults = sam2bam_singularity(mapResults)
    
    sortedBam = sortBam_singularity(bamResults)
    
    indexedBam = indexBam_singularity(sortedBam)
    
    bedGr = makeGenomeCov_singularity(sortedBam, reference_path)

    bedGr.view { "Coverage Plot: ${it}" }

    callPerl(bedGr).view()   
}
