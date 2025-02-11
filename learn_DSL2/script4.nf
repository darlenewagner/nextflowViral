#! /apps/x86_64/nextflow/23.10.0 nextflow

nextflow.enable.dsl=2

// Counts reads in paired fastq.gz file

params.inputPair = "${baseDir}/bowtieConsensTestFiles/eng_live_atten_poliovirus/polio_sample_3_screened_trim_R?_001.fastq.gz"
params.outdir = "${baseDir}/learn_DSL2/local_output/"

process READCOUNT {

    publishDir "${baseDir}/learn_DSL2/local_output/", mode: 'copy'

    input:
    tuple val(sample_name), path(reads)

    output:
    path "${sample_name}.readcounts.txt"

    script:
    """
    echo ${sample_name} > ${sample_name}.readcounts.txt
    zcat ${reads[0]} | awk 'END {print NR/4}' >> ${sample_name}.readcounts.txt
    zcat ${reads[1]} | awk 'END {print NR/4}' >> ${sample_name}.readcounts.txt
    """
}

workflow {
    Channel
        .fromFilePairs(params.inputPair, checkIfExists: true)
        .set { read_pairs_ch }

    READCOUNT(read_pairs_ch)

}
