# nextflowViral
All pipelines are written in Nextflow and/or Python to facilitate comparative genomics analysis for viruses.  All these pipelines require Nextflow 21.04.3 and above and most of the pipelines require Python 3.9.1 and above.

## bcfConsensusFromPairedFastq.nf
A Nextflow pipeline for mapping an arbitrary number of paired-end files to a viral reference genome to produce consensus sequences in fasta format

```nextflow run bcfConsensusFromPairedFastq.nf```

Prerequisites: smalt, samtools, and bcftools
Input folders: paired_files/ and Calicivirus_references/
Output folders: mapOut/ and multi_consensus/
