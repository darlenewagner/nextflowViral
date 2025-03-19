#! /apps/x86_64/nextflow/24.04.2 nextflow

nextflow.enable.dsl=2

// Counts reads in paired fastq.gz file

params.inputPair = "${baseDir}/../../bowtieConsensTestFiles/eng_live_atten_poliovirus/polio-sample-8_S13_R{1,2}_001.fastq.gz"
params.outdir = "${baseDir}/../../output/"

process READCOUNT {

    publishDir "${baseDir}/../../output/", mode: 'copy'

    input:
    tuple val(sample_name), path(reads)

    output:
    tuple val(sample_name), path("${sample_name}.readcounts.txt")

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

    see = READCOUNT(read_pairs_ch)

    see.view {"Read counts: ${it}"}

}
