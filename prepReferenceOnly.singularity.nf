#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/* Usage
     nextflow run prepReferenceOnly.singularity.nf --reference <fasta file>
     Web-connection-dependent Singularity containerization which calls bowtie2 and samtools
 */

params.reference = "${baseDir}/bowtieConsensTestFiles/eng_live_atten_poliovirus/MZ245455.1"



process bowtie2map_singularity {

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/bowtie2:2.5.4--he96a11b_5' :
    'quay.io/biocontainers/bowtie2:2.5.4--he96a11b_5' }"
    
    //publishDir "${baseDir}/bowtieConsensTestFiles/eng_live_atten_poliovirus/", mode: 'copy'

    input:
    path reference
    
    output:
    tuple val(reference), path("${reference}.*")
    
    script:
    """
    bowtie2-build 
    """
}

    