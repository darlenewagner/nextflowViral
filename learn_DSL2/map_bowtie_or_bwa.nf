#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/* Usage
     nextflow run word_count.nf --inputPair <input_file>
 */

params.inputPair = "data/untrimmed_fastq/SRR2584863_1.fastq.gz"
