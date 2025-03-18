# nextflowViral
This set of pipelines represents a viral genomic SNPs and consensus computing utility with a minimum number of software dependencies compared to comparable tools.

### Running local HPC Cluster consensusFromPairedFastq.nf

```module load nextflow/24.04.2````

```module load bowtie2/2.3.5.1```

```module load samtools/1.9```

```module load bcftools/1.9```

```module load htslib/1.19.1```

```module load bedtools/2.27.1```

```module load seqtk/1.3```


```nextflow run consensHiCovFromPairedFastq.nf --reference "$PWD/bowtieConsensTestFiles/EnterovirusD70/MT081369_JPN_1989-23292" --inputPair "$PWD/bowtieConsensTestFiles/EnterovirusD70/EnterovirusD70_SRR13402413_R{1,2}_001.fastq.gz"```

After successful run, folder *output/* should contain *EnterovirusD70_SRR13402413_Rb.bedGraph*, *\*.bedGraph.FiveX*, *\*.fasta*, *\*.fiveX.fasta*, *\*.snp.tsv*, *\*.vcf*, *\*.vcf.gz*, and **\.vcf.gz.csi*. 

To run without parameters, be prepared to update parameters, inputPair and reference, in ***nextflow.config***.



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
