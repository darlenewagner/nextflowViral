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
    tuple val(sample_name), path(reads)
    // val(reference)

    output:
    stdout

    script:
    """
     if [[ "${reads[0]}" =~ .*gz\$ ]];
     then
        printf '${reads[0]} first line is: ->  '
        zcat ${reads[0]} | head -1
     else
        printf '${reads[0]} first line is ->  '
        head -1 ${reads[0]}
     fi
    """
  
  }


workflow {
   // Channel
   //     .fromFilePairs(params.inputPair, checkIfExists: true)
   //     .set { read_pairs_ch }
   // Channel.fromFile(params.reference, checkIfExists: true)
   //       .set { ref }

    Channel
        .fromFilePairs(params.inputPair, checkIfExists: true)
        .set { read_pairs_ch }


    LOOKSY(read_pairs_ch) | view
    
   // LOOKSY.out.view()    

    // bowtie2map(read_pairs_ch, params.ref)

}