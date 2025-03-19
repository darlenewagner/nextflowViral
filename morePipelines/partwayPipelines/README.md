## What partwayPipeline scripts do and how to run them
### - The pipelines in this subfolder were written to perform simple, 1-step or 2-step tasks
#### *Installing Prerequisites:*
      In parent directory:
      ```module load python/3.10.8```
      ```python Python_scripts/prepReferenceOnly.py bowtieConsensTestFiles/eng_live_atten_poliovirus/```
      ```cd morePipelines/partwayPipelines/```
### - script4.nf counts reads in paired end .fastq.gz files

### - map_bowtie_or_bwa.nf maps test set reads, polio-sample-8_S13_R{1,2}_001.fastq.gz to MZ245455
