# nextflowViral
All pipelines are written in Nextflow and/or Python to facilitate comparative genomics analysis for viruses.  All these pipelines require Nextflow 21.04.3 and above and most of the pipelines require Python 3.9.1 and above.

### bcfConsensusFromPairedFastq.nf
A Nextflow pipeline for mapping an arbitrary number of paired-end files to a viral reference genome to produce consensus sequences in fasta format

```nextflow run bcfConsensusFromPairedFastq.nf```

**Prerequisites:** smalt, samtools, and bcftools<br>
*Required input folders* are paired_files/ and Calicivirus_references/<br>
*Expected output folders* are mapOut/ and multi_consensus/

### bowtieConsensusFromPairedFastq.nf
Maps an arbitrary number of paired-end files to a viral reference genome to produce consensus sequences in fasta format and genome coverage files in .bed format

```nextflow run bowtieConsensusFromPairedFastq.nf```

**Prerequisites:** bowtie2, samtools, bedtools, and bcftools<br>
*Required input folders* are paired_files/ and ref_180330_Calici36_SureSelectNoro_Capture/<br>
***ref_180330_Calici36_SureSelectNoro_Capture*** contains .fasta, .sizes, and bowtie2-build .bt2 files
*Expected output folders* are bowtieOut/ and new_multi_consensus/
