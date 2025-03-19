## How to map test single-end reads to test reference
### Vignette 1: HPC Cluster installation in nextflowViral folder

`module load python/3.12.3`

`module load bowtie2/2.3.5.1`

`module load samtools/1.9`

`module load bedtools/2.27.1`

`module load nextflow/24.04.2`

#### Note: Assume tar -xfzv bowtieConsensTestFiles.tar.gz has already been run

`python Python_scripts/prepReferenceOnly.py bowtieConsensTestFiles/adenovirus_B3/`

`cd learn_DSL2/`

`nextflow run mapSingle_bowtie_to_plot_coverage.nf --reference "$PWD/../bowtieConsensTestFiles/adenovirus_B3/OR777202.1" --inputSingle "$PWD/../bowtieConsensTestFiles/adenovirus_B3/Pool-1_S1_adenovirus_B3_001.fastq.gz" --intermediate "$PWD/intermediate/"`

---

### Coming Soon: Vignette 2, Singularity containerized installation

