#! /apps/x86_64/nextflow/23.10.0 nextflow

nextflow.enable.dsl=2

// Counts reads in unpaired fastq file

params.inputSingle = "${baseDir}/../../bowtieConsensTestFiles/adenovirus_B3/Pool-1_S1_adenovirus_B3_001.fastq.gz"
params.outdir = "${baseDir}/../../output/"

process READCOUNT {

    publishDir "${baseDir}/../../output/", mode: 'copy'

    //tag "$sample.name"

    input:
    path sample_name

    output:
    tuple val(sample_name), path("${sample_name.simpleName}.readcounts.txt")
    
    script:
    """
    echo ${sample_name} > ${sample_name.simpleName}.readcounts.txt
    zcat ${sample_name} | awk 'END {print NR/4}' >> ${sample_name.simpleName}.readcounts.txt
    """
}

workflow {
    Channel
        .fromPath(params.inputSingle, checkIfExists: true)
        .set { fastq_ch }

   see = READCOUNT(fastq_ch)

   see.view { "Read counts: {$it}" }
   
}
