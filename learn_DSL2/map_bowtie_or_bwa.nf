#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/* Usage
     nextflow run word_count.nf --reference <fasta file> --inputPair <paired fastq files>
 */

params.reference = "${baseDir}/bowtieConsensTestFiles/eng_live_atten_poliovirus/MZ245455.1.fasta"

params.inputPair = "${baseDir}/bowtieConsensTestFiles/eng_live_atten_poliovirus/polio_sample_3_screened_trim_R?_001.fastq.gz"

process LOOKSY  // for debugging and sanity checking
  {

    input:
    tuple path(seeFile)

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
   // Channel
   //     .fromFilePairs(params.inputPair, checkIfExists: true)
   //     .set { read_pairs_ch }
   // Channel.fromFile(params.reference, checkIfExists: true)
    
    LOOKSY(params.inputPair) | view
    
   // LOOKSY.out.view()    

    // bowtie2map(read_pairs_ch, params.ref)

}