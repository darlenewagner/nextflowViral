#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/* Usage
     nextflow run prepReferenceOnly.singularity.nf --reference <fasta file>
     Web-connection-dependent Singularity containerization which calls bowtie2 and samtools

params.reference = "${baseDir}/bowtieConsensTestFiles/eng_live_atten_poliovirus/MZ245455.1.fasta"
params.outdir = "${baseDir}/bowtieConsensTestFiles/eng_live_atten_poliovirus/"
reference_name = file(params.reference).name
reference_idx = "${reference_name}.fai"
reference_path = file(params.reference).parent
*/

params.reference = "${baseDir}/bowtieConsensTestFiles/eng_live_atten_poliovirus/MZ245455.1.fasta"
reference_path = file(params.reference).parent

process buildBowtie2Index {

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/bowtie2:2.5.4--he96a11b_5' :
    'quay.io/biocontainers/bowtie2:2.5.4--he96a11b_5' }"

    publishDir "${baseDir}/bowtieConsensTestFiles/eng_live_atten_poliovirus/", mode: 'copy'

    input:
    tuple val(fasta), path(fasta_file)

    //output:
    //tuple val(fasta), path("${fasta_file}.*bt2")

    """
    bowtie2-build "${reference_path}"/"${fasta_file}.fasta" "${reference_path}"/"${fasta}"
    """
}

workflow {

   // params.fasta = file(params.fasta ?: 'bowtieConsensTestFiles/eng_live_atten_poliovirus/MZ245455.1.fasta')

   Channel.fromPath( params.reference )
          .map { file -> tuple(file.baseName, file) }
	  .set { fasta_file }

     buildBowtie2Index( fasta_file )
 }
