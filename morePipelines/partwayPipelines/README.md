## What partwayPipeline scripts do and how to run them
### - The pipelines in this subfolder were written to perform simple, 1-step or 2-step tasks
#### *Installing Prerequisites:*

In parent directory:

```module load python/3.10.8```

```python Python_scripts/prepReferenceOnly.py bowtieConsensTestFiles/eng_live_atten_poliovirus/```

```cd morePipelines/partwayPipelines/```

### - script4.nf counts reads in paired end .fastq.gz files

```module load nextflow/24.04.2```

```nextflow run script4.nf```

To see the result, look in *../../output/* for the file, *polio-sample-8_S13_R.readcounts.txt*.

To change inputPairs parameter, edit *nextflow.config* in the current folder, *morePipelines/partwayPipelines/*.

Alternatively, provide an input path as an argument to *script4.nf* as shown:

```nextflow run script4.nf --inputPair "$PWD/../../myViralReads/virus_sample_X_R{1,2}_001.fastq.gz"```

### - map_bowtie_or_bwa.nf maps reads to an indexed reference 
#### *Test reads polio-sample-8_S13_R{1,2}_001.fastq.gz to test reference MZ245455*

```module load bowtie2/2.4.5```

```nextflow run map_bowtie_or_bwa.nf```

#### *User-defined reads and reference*

```python ../../Python_scripts/prepReferenceOnly.py ../../bowtieConsensTestFiles/EnterovirusD70/```

```nextflow run map_bowtie_or_bwa.nf --reference "$PWD/../../bowtieConsensTestFiles/EnterovirusD70/MT081369_JPN_1989-23292" --inputPair "$PWD/../../bowtieConsensTestFiles/EnterovirusD70/EnterovirusD70_SRR13402413_R{1,2}_001.fastq.gz"```

