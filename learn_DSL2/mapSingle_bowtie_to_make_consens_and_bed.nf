#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/* Usage
     nextflow run word_count.nf --reference <fasta file> --inputPair <paired fastq files>
 */

params.reference = "${baseDir}/bowtieConsensTestFiles/eng_live_atten_poliovirus/MZ245455.1"
reference_name = file(params.reference).name
reference_idx = "${reference_name}.fai"
reference_path = file(params.reference).parent

params.inputSingle = "${baseDir}/../bowtieConsensTestFiles/pseudoPairs/Pool-1_S1_final_pseudoPairs.fastq"
inputSingle_path = file(params.inputSingle).parent
params.output = "${baseDir}/../output/"

process LOOKSY  // for debugging and sanity checking
  {

    input:
    path(sample_name)
    path(reference)
    
    output:
    stdout
    
    script:
    """
     printf '${reference} first line is: -> '
     head -1 ${reference}
     if [[ "${sample_name}" +~ .*fastq\$ ]];
     then
        printf '${sample_name} first line is: ->  '
        zcat ${sample_name} | head -1
        printf '${sample_name} first line is: ->  '
        cat ${sample_name} | head -1
     fi
    """
  
  }


process bowtie2map {

    publishDir "${baseDir}/../intermediate/", mode: 'copy'
    
    input:
    path(sample_name)
    path reference 

    output:
    tuple val(sample_name), path("${sample_name}.sam")    

    script:
    """
    bowtie2 --no-unal --no-mixed -x "${reference}"/"${reference_name}" -U "${sample_name}" > "${sample_name}.sam"
    """
}


process bowtie2map_singularity {

    tag { sample }
    
    publishDir "${baseDir}/../intermediate/", mode: 'copy'
    
    input:
    path(sample_name)
    path reference 
    
    output:
    tuple val(sample_name.baseName), path("${sample_name.baseName}.sam")    
       
    script:
    """
    singularity exec "${baseDir}/../"my_bowtie2.sif bowtie2 --no-unal --no-mixed -x "${reference}"/"${reference_name}" -U "${sample_name}" > "${sample_name.baseName}.sam"
    """
}


process sam2bam {

    publishDir "${baseDir}/../intermediate/", mode: 'copy'
    
    input:
    tuple val(sample_name), path("${sample_name}.sam")
    
    output:
    tuple val(sample_name), path("${sample_name}.bam")

    script:
    """
    samtools view "${sample_name}".sam -o "${sample_name}".bam
    """

}

process sortBam {

    publishDir "${baseDir}/../intermediate/", mode: 'copy'

    input:
    tuple val(sample_name), path("${sample_name}.bam")
    
    output:
    tuple val(sample_name), path("${sample_name}.sorted.bam")
    
    script:
    """
    samtools sort "${sample_name}".bam -o "${sample_name}".sorted.bam
    """
}

process indexBam {
    
    publishDir "${baseDir}/../intermediate/", mode: 'copy'
    
    input:
    tuple val(sample_name), path("${sample_name}.sorted.bam")
    
    output:
    tuple val(sample_name), path("${sample_name}.sorted.bam.bai")
    
    script:
    """
    samtools index -b "${sample_name}".sorted.bam "${sample_name}".sorted.bam.bai
    """

}

process makeVCF {

    publishDir "${baseDir}/../output/", mode: 'copy'

    input:
    tuple val(sample_name), path("${sample_name}.sorted.bam")
    path reference
    
    output:
    tuple val(sample_name), path("${sample_name}.vcf")
    
    script:
    """
    bcftools mpileup -d 35000 -Ob -f "${reference}"/"${reference_name}".fasta -Q 20 -q 20 --annotate FORMAT/AD,FORMAT/DP "${sample_name}".sorted.bam | bcftools call -mv --ploidy 1 -Ov --output "${sample_name}".vcf
    """
    
}

process zipVCF {
    
    publishDir "${baseDir}/../output/", mode: 'copy'
    
    input:
    tuple val(sample_name), path("${sample_name}.vcf")
    
    output:
    tuple val(sample_name), path("${sample_name}.vcf.gz")

    
    script:
    """
      bgzip -c "${sample_name}".vcf > "${sample_name}".vcf.gz  

    """
}

process csiVCF {

   publishDir "${baseDir}/../output/", mode: 'copy'

    input:
    tuple val(sample_name), path("${sample_name}.vcf.gz")

    output:
    tuple val(sample_name), path("${sample_name}.vcf.gz.csi")

    script:
    """
      bcftools index "${sample_name}".vcf.gz -o "${sample_name}".vcf.gz.csi    
    """
 }

process makeBcfConsensus {
    
    publishDir "${baseDir}/../output/", mode: 'copy'
    
    input:
    tuple val(sample_name), path("${sample_name}.vcf.gz")
    tuple val(sample_name), path("${sample_name}.vcf.gz.csi")
    path reference
    
    output:
    tuple val(sample_name), path("${sample_name}.fasta")
    
    script:
    """
    cat "${reference}"/"${reference_name}".fasta | bcftools consensus "${sample_name}".vcf.gz --sample "${sample_name}".sorted.bam -o "${sample_name}".fasta
    """
    
}


workflow {
    
    Channel
        .fromPath(params.inputSingle, checkIfExists: true)
        .set { fastq_ch }
   
   
   // LOOKSY(fastq_ch, params.reference) | view
    
   //  mapResults = bowtie2map(fastq_ch, reference_path)
    
    mapResults = bowtie2map_singularity(params.inputSingle, reference_path) 
    
    mapResults.view { "Bowtie2 Results: ${it}" }
    
    bamResults = sam2bam(mapResults)
    
    sortedBam = sortBam(bamResults)
    
    indexedBam = indexBam(sortedBam)
    
    myVCF = makeVCF(sortedBam, reference_path)   
    
    myVCF.view { "SNP calls: ${it}" }
    
    myVCFz = zipVCF(myVCF)

    myVCFcsi = csiVCF(myVCFz)

    fasta = makeBcfConsensus(myVCFz, myVCFcsi, reference_path)

    fasta.view { "Consensus FASTA: ${it}" }
}
