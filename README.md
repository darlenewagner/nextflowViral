# nextflowViral
This set of pipelines represents a viral genomic SNPs and consensus computing utility with a minimum number of software dependencies compared to comparable tools.

### Vignette 1: Install and run from LINUX/UNIX using Singularity and Nextflow only - *Recommended Method*



### Vignette 2: Install prereqs for running local HPC Cluster consensusFromPairedFastq.nf

```module load nextflow/24.04.2```

```module load bowtie2/2.3.5.1```

```module load samtools/1.9```

```module load bcftools/1.9```

```module load htslib/1.19.1```

```module load bedtools/2.27.1```

```module load seqtk/1.3```

#### Vignette 2a: No command line arguments. Input files defined in *nextflow.config*

```nextflow run consensHiCovFromPairedFastq.nf```

After successful run, folder *output/* should contain *polio-sample-8_S13_R.bedGraph*, *\*.bedGraph.FiveX*, *\*.fasta*, *\*.fiveX.fasta*, *\*.snp.tsv*, *\*.vcf*, *\*.vcf.gz*, and **\.vcf.gz.csi*. The file, *polio-sample-8_S13_R.snp.tsv*, should contain 1 SNP position. 

#### Vignette 2b: Alternate --reference and --inputPair defined as input arguments.

```nextflow run consensHiCovFromPairedFastq.nf --reference "$PWD/bowtieConsensTestFiles/EnterovirusD70/MT081369_JPN_1989-23292" --inputPair "$PWD/bowtieConsensTestFiles/EnterovirusD70/EnterovirusD70_SRR13402413_R{1,2}_001.fastq.gz"```

After successful run, folder *output/* should contain *EnterovirusD70_SRR13402413_R.bedGraph*, *\*.bedGraph.FiveX*, *\*.fasta*, *\*.fiveX.fasta*, *\*.snp.tsv*, *\*.vcf*, *\*.vcf.gz*, and **\.vcf.gz.csi*. The file, *EnterovirusD70_SRR13402413_R.snp.tsv*, should not contain SNP positions. 

To point parameters to other files, be prepared to update parameters, inputPair and reference, in ***nextflow.config***.




