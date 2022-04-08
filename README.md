# nextflowViral
All pipelines are written in Nextflow and/or Python to facilitate comparative genomics analysis for viruses.  All these pipelines require Nextflow 21.04.3 and above and most of the pipelines require Python 3.9.1 and above.

### bcfConsensusFromPairedFastq.nf
A Nextflow pipeline for mapping an arbitrary number of paired-end files to a viral reference genome to produce consensus sequences in fasta format

```nextflow run bcfConsensusFromPairedFastq.nf```

<u>Prerequisites:</u> smalt, samtools, and bcftools
Required input folders: paired_files/ and Calicivirus_references/
Expected output folders: mapOut/ and multi_consensus/

### bowtieConsensusFromPairedFastq.nf
Maps an arbitrary number of paired-end files to a viral reference genome to produce consensus sequences in fasta format and genome coverage files in .bed format

```nextflow run bowtieConsensusFromPairedFastq.nf```

Required input folders: paired_files/ and Calicivirus_references/
Expected output folders: bowtieOut/ and new_multi_consensus/
