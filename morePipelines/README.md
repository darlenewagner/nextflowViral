## Vignette 1: How to map paired-end reads to reference using locally-installed prerequisites

`module load nextflow/24.04.2`

`module load bowtie2/2.3.5.1`

`module load samtools/1.9`

`module load bcftools/1.9`

`module load htslib/1.19.1`

`module load bedtools/2.27.1`

`module load seqtk/1.3`

`nextflow run prepReferenceOnly.singularity.nf --makeReference "$PWD/bowtieConsensTestFiles/eng_live_atten_poliovirus/MZ245455.1.fasta" --indexdir "$PWD/bowtieConsensTestFiles/eng_live_atten_poliovirus/"`

`cd morePipelines/`

`nextflow run consensHiCovFromPairedFastq.nf`

#### Note: Presence of programs required by processes result in printing 'true' or 'false' accordingly

## Vignette 2: How to map single-end reads to reference for breadth and depth of coverage

### Vignette 2a: Singularity containerized installation and run

#### In *nextflowViral/* home folder:

`nextflow run prepReferenceOnly.singularity.nf --makeReference "$PWD/bowtieConsensTestFiles/adenovirus_B3/OR777202.1.fasta" --indexdir "$PWD/bowtieConsensTestFiles/adenovirus_B3/"`

`cd morePipelines/`

#### In *morePipelines/*

`nextflow run mapSingle_bowtie_to_plot_coverage.singularity.nf --reference "$PWD/../bowtieConsensTestFiles/adenovirus_B3/OR777202.1" --inputSingle "$PWD/../bowtieConsensTestFiles/adenovirus_B3/Pool-1_S1_adenovirus_B3_001.fastq.gz" --intermediate "$PWD/../intermediate/" --output "$PWD/../outdir/"`

Final line of output to screen:
*OR777202.1      Avg. Depth: 111.8       Max. Breadth: 81.9%*


### Vignette 2b: HPC Cluster installation in nextflowViral folder

`module load python/3.12.3`

`module load bowtie2/2.3.5.1`

`module load samtools/1.9`

`module load bedtools/2.27.1`

`module load nextflow/24.04.2`

#### Note: Assume tar -xfzv bowtieConsensTestFiles.tar.gz has already been run

`python Python_scripts/prepReferenceOnly.py bowtieConsensTestFiles/adenovirus_B3/`

`cd morePipelines/`

`nextflow run mapSingle_bowtie_to_plot_coverage.nf --reference "$PWD/../bowtieConsensTestFiles/adenovirus_B3/OR777202.1" --inputSingle "$PWD/../bowtieConsensTestFiles/adenovirus_B3/Pool-1_S1_adenovirus_B3_001.fastq.gz" --intermediate "$PWD/intermediate/"`



