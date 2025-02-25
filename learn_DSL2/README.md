### Running test pipeline, script4.nf

`module load nextflow/24.04.2`

`nextflow run learn_DSL2/script4.nf`

### Running the Nextflow and Bowtie test pipeline:

`module load nextflow/24.04.2`

`module load bowtie2/2.3.5.1`

`nextflow run learn_DSL2/map_bowtie_or_bwa.nf --reference "$PWD/bowtieConsensTestFiles/eng_live_atten_poliovirus/MZ245455.1" --inputPair "$PWD/bowtieConsensTestFiles/eng_live_atten_poliovirus/polio_sample_3_screened_trim_R{1,2}_001.fastq.gz"`

### Running the Nextflow, Bowtie, Samtools, and Bcftools test pipeline:

`module load nextflow/24.04.2`

`module load samtools/1.9`

`module load bcftools/1.9`

`nextflow run learn_DSL2/map_bowtie_to_vcf.nf --reference "$PWD/bowtieConsensTestFiles/eng_live_atten_poliovirus/MZ245455.1"`

### Coming Soon: Singularity Bowtie2

`module load nextflow/24.04.2`

`python Singularity/get_bowtie2.py`

`nextflow run learn_DSL2/map_bowtie_or_bwa.nf --reference "$PWD/bowtieConsensTestFiles/eng_live_atten_poliovirus/MZ245455.1"`