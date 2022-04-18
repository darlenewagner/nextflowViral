# nextflowViral
The pipeline, *bowtieConsensusFromPairedFastq.nf*, generates consensus genomes from paired fastq files mapped to reference at minimum 1X and 25X coverage.  To run, it requires Nextflow 21.04.3 and above and for setup, it requires Python 3.9.1 and above.

### Running bowtieConsensusFromPairedFastq.nf
To run using input parameters in *nextflow.config*, simply type the following after installation:

```nextflow run bowtieConsensusFromPairedFastq.nf```

To run using custom parameters for *--querydir* and *--reference*, for example, the folder, my_paired_fastq/, and the reference basename, my_reference/my_genome, type:

**Prerequisites:** bowtie2, samtools, bedtools, and bcftools<br>
*Required input folders* are paired_files/ and ref_180330_Calici36_SureSelectNoro_Capture/<br>
***ref_180330_Calici36_SureSelectNoro_Capture/*** contains .fasta, .sizes, and bowtie2-build .bt2 files<br>
*Expected output folders* are bowtieOut/ and new_multi_consensus/

### Installation and setup for bowtieConsensusFromPairedFastq.nf ###
