# nextflowViral
In the nextflowViral package, the pipeline script, ***bowtieConsensHiCovFromPairedFastq.nf***, is written in Nextflow DSL 1 and is thus out-of-date. It generates consensus genomes from paired fastq files mapped to reference at minimum 1X and 25X coverage.  To run old_DSL1 files, Nextflow 21.04.3 is required.  To run the new pipelines compatible with DSL 2 (under learn_DSL2), Nextflow 23.10.0 and above is required.  For setup, Python 3.9.1 and above is required.

### Running bowtieConsensusFromPairedFastq.nf
To run using input parameters in ***nextflow.config***, simply type the following after installation:

```nextflow run bowtieConsensHiCovFromPairedFastq.nf```

To run using custom parameters for *--querydir* and *--reference*, for example, the folder, *my_paired_fastq/*, and the reference basename, *my_reference/my_genome*, type:

```nextflow run bowtieConsensHiCovFromPairedFastq.nf --querydir "$PWD/my_paired_fastq/*_R{1,2}*.fastq" --reference $PWD/my_reference/my_genome```

Alternatively, to avoid typing complicated filepaths at the command line, simply edit the *querydir*, *reference*, *intermdir*, and/or *consensdir* in the ***nextflow.config*** file.

**Prerequisites:** bowtie2.2.3.X.X, samtools/1.15.X, BEDTools/2.26.X, bcftools/1.10.X, seqtk/1.3, perl/5.16.X-MT, and nextflow/21.04.X<br>
*Required/default input folders* are D70_paired_files/ and D70_reference_genome/<br>
***ref_180330_Calici36_SureSelectNoro_Capture/*** contains .fasta, .sizes, and bowtie2-build .bt2 files<br>
*Expected/default output folders* are bowtieConsensInterm/ and bowtieConsensOutput/

### Installation and setup for bowtieConsensusFromPairedFastq.nf ###
To set up run environment with required input and output folders and appropriately-formatted files within input folders, enter the command:<br><br>
```python setupFolders.py --input newTest.tar```
<br><br>
Here, *--input newTest.tar* is optional.  The default input is *bowtieConsensTestFiles.tar*.
