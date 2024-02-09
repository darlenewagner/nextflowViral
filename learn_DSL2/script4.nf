#! /apps/x86_64/nextflow/23.10.0 nextflow

nextflow.enable.dsl=2

params.reads = "$projectDir/Polio_MiSeq_with_QC/FirstFive_vs_3015815424/*_R{1,2}_001.fastq"
params.reference = "$projectDir/NCBI_Polio_reference/AY184220.1"
params.outdir = "$projectDir/redevelop_Interm_2024"

process READCOUNT {

//    tag { see }
    publishDir params.outdir, mode:"copy"

    input:
    tuple val(sample_name), path(reads)

    output:
    path "${sample_name}.txt"

    script:
    """
    echo ${sample_name} > ${sample_name}.txt
    awk 'END {print NR/4}' ${reads[0]} >> ${sample_name}.txt
    awk 'END {print NR/4}' ${reads[1]} >> ${sample_name}.txt
    """
}

workflow {
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs_ch }

    READCOUNT(read_pairs_ch)

}
