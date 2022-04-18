# nextflowViral
The pipeline, *bowtieConsensusFromPairedFastq.nf*, generates consensus genomes from paired fastq files mapped to reference at minimum 1X and 25X coverage.  To run, it requires Nextflow 21.04.3 and above and for setup, it requires Python 3.9.1 and above.

### Running bowtieConsensusFromPairedFastq.nf
Maps an arbitrary number of paired-end files to a viral reference genome to produce consensus sequences in fasta format and genome coverage files in .bed format

```nextflow run bowtieConsensusFromPairedFastq.nf```

**Prerequisites:** bowtie2, samtools, bedtools, and bcftools<br>
*Required input folders* are paired_files/ and ref_180330_Calici36_SureSelectNoro_Capture/<br>
***ref_180330_Calici36_SureSelectNoro_Capture/*** contains .fasta, .sizes, and bowtie2-build .bt2 files<br>
*Expected output folders* are bowtieOut/ and new_multi_consensus/

### Setup for bowtieConsensusFromPairedFastq.nf ###
