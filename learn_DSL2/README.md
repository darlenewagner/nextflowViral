### Running the Nextflow and Bowtie test pipeline:

`module load nextflow/24.04.2`

`module load bowtie2/2.3.5.1`

`nextflow run learn_DSL2/map_bowtie_or_bwa.nf --reference "$PWD/bowtieConsensTestFiles/eng_live_atten_poliovirus/MZ245455.1" --inputPair "$PWD/bowtieConsensTestFiles/eng_live_atten_poliovirus/polio_sample_3_screened_trim_R{1,2}_001.fastq.gz"`

### Coming Soon: Singularity Bowtie2
