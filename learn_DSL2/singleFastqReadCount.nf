#! /apps/x86_64/nextflow/23.10.0 nextflow

nextflow.enable.dsl=2

// Counts reads in unpaired fastq file

params.inputSingle = "${baseDir}/../bowtieConsensTestFiles/pseudoPairs/Pool-1_S1_final_pseudoPairs.fastq"
params.outdir = "${baseDir}/learn_DSL2/local_output/"

process READCOUNT {

    publishDir "${baseDir}/learn_DSL2/local_output/", mode: 'copy'

    input:
    path(sample_name)

    output:
    path "${sample_name}.readcounts.txt"
    
    
    script:
    """
    echo ${sample_name} > ${sample_name}.readcounts.txt
    cat ${sample_name} | awk 'END {print NR/4}' >> ${sample_name}.readcounts.txt
    """
}

workflow {
    Channel
        .fromPath(params.inputSingle, checkIfExists: true)
        .set { fastq_ch }

    READCOUNT(fastq_ch)

}
