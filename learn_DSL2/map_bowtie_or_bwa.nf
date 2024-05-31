#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/* Usage
     nextflow run word_count.nf --reference <fasta file> --inputPair <paired fastq files>
 */

params.reference = "$PWD/../bowtieConsensTestFiles/Sabin-2_reference/AY184220.1.fasta"

params.inputPair = "$PWD/../../nextflow_2023_for_read_mapping/Polio_MiSeq_trimmomatic/LastEight_vs_3015821190/3015821227_S13_R?_001.fastq.gz"

process LOOKSY  // for debugging and sanity checking
  {

    input:
    path seeFile

    output:
    stdout

    script:
    """
     if [[ "${seeFile}" =~ .*gz\$ ]];
     then
        printf '${seeFile} first line is: ->  '
        gunzip -c ${seeFile} | head -1
     else
        printf '${seeFile} first line is ->  '
        head -1 ${seeFile}
     fi
    """
  
  }


workflow {
    Channel
        .fromFilePairs(params.inputPair, checkIfExists: true)
        .set { read_pairs_ch }

    Channel.fromFile(params.reference, checkIfExists: true)
    
    LOOKSY(params.reference)
    
    LOOKSY.out.view()    

    // bowtie2map(read_pairs_ch, params.ref)

}