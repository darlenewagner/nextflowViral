#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/* Usage
     nextflow run prepReferenceOnly.singularity.nf --makeReference <fasta file> --indexDir <destination folder>
     Web-connection-dependent Singularity containerization which calls bowtie2 and samtools
*/

params.makeReference = "${baseDir}/bowtieConsensTestFiles/eng_live_atten_poliovirus/MZ245455.1.fasta"
params.indexDir = "${baseDir}/bowtieConsensTestFiles/eng_live_atten_poliovirus/"
reference_path = file(params.makeReference).parent

process buildBowtie2Index {

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/bowtie2:2.5.4--he96a11b_5' :
    'quay.io/biocontainers/bowtie2:2.5.4--he96a11b_5' }"

    publishDir "${params.indexDir}", mode: 'copy'

    input:
    tuple val(fasta), path(fasta_file)

    """
    bowtie2-build "${reference_path}"/"${fasta_file}" "${reference_path}"/"${fasta}"
    """
}

process samtoolsFaidx {

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/samtools:1.9--h91753b0_8' :
    'quay.io/biocontainers/samtools:1.9--h91753b0_8' }"

    
    publishDir "${reference_path}", mode: 'copy'

    input:
    tuple val(fasta), path(fasta_file)

    output:
    stdout
   // tuple val(fasta), path("${fasta_file}.fai")    
   // tuple val(fasta), path("${fasta_file}.sizes")    

    """
    samtools faidx "${reference_path}"/"${fasta_file}"
    cut -f 1,2 "${reference_path}"/"${fasta_file}.fai"
    """

}

workflow {


   Channel.fromPath( params.makeReference )
          .map { file -> tuple(file.baseName, file) }
	  .set { fasta_file }

     buildBowtie2Index( fasta_file )
     
     samtoolsFaidx( fasta_file ).collectFile( name: "${fasta_file}.sizes" ).view()
     
  }
