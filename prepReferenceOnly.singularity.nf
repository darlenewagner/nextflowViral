#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/* Usage
     nextflow run prepReferenceOnly.singularity.nf --reference <fasta file>
     Web-connection-dependent Singularity containerization which calls bowtie2 and samtools
 */

params.reference = "${baseDir}/bowtieConsensTestFiles/eng_live_atten_poliovirus/MZ245455.1.fasta"
params.outdir = "${baseDir}/bowtieConsensTestFiles/eng_live_atten_poliovirus/"
reference_name = file(params.reference).name
reference_idx = "${reference_name}.fai"
reference_path = file(params.reference).parent


process bowtie2map_singularity {

    tag "{reference_name}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/bowtie2:2.5.4--he96a11b_5' :
    'quay.io/biocontainers/bowtie2:2.5.4--he96a11b_5' }"
    
    //publishDir params.outdir, mode: 'copy'

    input:
 //   path fasta_file, type: 'file'
    val reference_name


    output:
    path "${reference_name}.bt2", type: 'file'
    
    script:
    """
    bowtie2-build ${params.reference} ${reference_name}
    """
}

workflow {
  
 bowtie2map_singularity(fasta_file: params.reference, sample_id: reference_name)
  
}
