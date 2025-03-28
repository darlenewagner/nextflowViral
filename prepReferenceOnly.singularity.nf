#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/* Usage
     nextflow run prepReferenceOnly.singularity.nf --reference <fasta file>
     Web-connection-dependent Singularity containerization which calls bowtie2 and samtools
*/

params.makeReference = "${baseDir}/bowtieConsensTestFiles/eng_live_atten_poliovirus/MZ245455.1"
params.indexdir = "${baseDir}/bowtieConsensTestFiles/eng_live_atten_poliovirus/"
reference_path = file(params.makeReference).parent

process buildBowtie2Index {

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/bowtie2:2.5.4--he96a11b_5' :
    'quay.io/biocontainers/bowtie2:2.5.4--he96a11b_5' }"

    publishDir "${params.indexdir}", mode: 'copy'

    input:
    tuple val(fasta), path(fasta_file)

    """
    bowtie2-build "${reference_path}"/"${fasta_file}.fasta" "${reference_path}"/"${fasta}"
    """
}

workflow {


   Channel.fromPath( params.makeReference )
          .map { file -> tuple(file.baseName, file) }
	  .set { fasta_file }

     buildBowtie2Index( fasta_file )

  }
