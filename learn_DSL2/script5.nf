#! /apps/x86_64/nextflow/23.10.0 nextflow

nextflow.enable.dsl=2

params.reads = "$projectDir/Polio_MiSeq_with_QC/FirstFive_vs_3015815424/*_R{1,2}_001.fastq"
params.ref = "$projectDir/NCBI_Polio_reference/AY184220.1"
params.outdir = "$projectDir/redevelop_Interm_2024"

process bowtie2map {

    publishDir params.outdir, mode:"copy"

    input:
    tuple val(sample_name), path(reads)
    val(reference)

    output:
    path "${sample_name}.sam"

    script:
    """
    /apps/x86_64/bowtie2/bowtie2-2.3.5.1/bowtie2-2.3.5.1-linux-x86_64/bowtie2 --no-unal --no-mixed -x ${reference} -1 ${reads[0]} -2 ${reads[1]} -S ${sample_name}.sam
    """
}

workflow {
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs_ch }

    bowtie2map(read_pairs_ch, params.ref)

}

