#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/* Usage
     nextflow run consensHiCovFromPairedFastq.singularity.nf --reference <fasta file> --inputPair <paired fastq files>
     
     Web-connection-dependent Singularity containerization which calls required container separately within each process definition
 */

params.reference = "${baseDir}/bowtieConsensTestFiles/eng_live_atten_poliovirus/MZ245455.1"
reference_name = file(params.reference).name
reference_idx = "${reference_name}.fai"
reference_path = file(params.reference).parent

params.inputPair = "${baseDir}/bowtieConsensTestFiles/eng_live_atten_poliovirus/polio-sample-2_S7_R{1,2}_001.fastq.gz"
params.output = "${baseDir}/output/"

process LOOKSY  // for debugging and sanity checking
  {

    input:
    tuple val(sample_name), path(reads)
    path(reference)

    output:
    stdout

    script:
    """
     printf '${reference} first line is: -> '
     head -1 ${reference}
     if [[ "${reads[0]}" =~ .*gz\$ ]] && [[ "${reads[1]}" =~ .*gz\$ ]];
     then
        printf '${reads[0]} first line is: ->  '
        zcat ${reads[0]} | head -1
        printf '${reads[1]} first line is: ->  '
        zcat ${reads[1]} | head -1
     else
        printf '${reads[0]} first line is ->  '
        head -1 ${reads[0]}
        printf '${reads[1]} first line is: ->  '
        head -1 ${reads[1]} 
     fi
    """
  
  }


process bowtie2map_singularity {

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/bowtie2:2.5.4--he96a11b_5' :
    'quay.io/biocontainers/bowtie2:2.5.4--he96a11b_5' }"
    
    publishDir "${baseDir}/intermediate/", mode: 'copy'
    
    input:
    tuple val(sample_name), path(reads)
    path reference 

    output:
    tuple val(sample_name), path("${sample_name}.sam")    

    script:
    """
    bowtie2 --no-unal --no-mixed -x "${reference}"/"${reference_name}" -1 "${reads[0]}" -2 "${reads[1]}" > "${sample_name}.sam"
    """
}


process sam2bam_singularity {

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/samtools:1.9--h91753b0_8' :
    'quay.io/biocontainers/samtools:1.9--h91753b0_8' }"

    publishDir "${baseDir}/intermediate/", mode: 'copy'
    
    input:
    tuple val(sample_name), path("${sample_name}.sam")
    
    output:
    tuple val(sample_name), path("${sample_name}.bam")

    script:
    """
    samtools view "${sample_name}".sam -o "${sample_name}".bam
    """

}

process sortBam_singularity {

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/samtools:1.9--h91753b0_8' :
    'quay.io/biocontainers/samtools:1.9--h91753b0_8' }"
    
    publishDir "${baseDir}/intermediate/", mode: 'copy'

    input:
    tuple val(sample_name), path("${sample_name}.bam")
    
    output:
    tuple val(sample_name), path("${sample_name}.sorted.bam")
    
    script:
    """
    samtools sort "${sample_name}".bam -o "${sample_name}".sorted.bam
    """
}

process indexBam_singularity {
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/samtools:1.9--h91753b0_8' :
    'quay.io/biocontainers/samtools:1.9--h91753b0_8' }"
    
    publishDir "${baseDir}/intermediate/", mode: 'copy'
    
    input:
    tuple val(sample_name), path("${sample_name}.sorted.bam")
    
    output:
    tuple val(sample_name), path("${sample_name}.sorted.bam.bai")
    
    script:
    """
    samtools index -b "${sample_name}".sorted.bam "${sample_name}".sorted.bam.bai
    """

}

process makeVCF_singularity {

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/bcftools:1.9--ha228f0b_4' :
    'quay.io/biocontainers/bcftools:1.9--ha228f0b_4' }"
    
    publishDir "${baseDir}/output/", mode: 'copy'

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

process zipVCF_singularity {
    
   container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
   'https://depot.galaxyproject.org/singularity/htslib:1.19.1--h81da01d_2' :
   'quay.io/biocontainers/htslib:1.19.1--h81da01d_2' }"
    
    publishDir "${baseDir}/output/", mode: 'copy'
    
    input:
    tuple val(sample_name), path("${sample_name}.vcf")
    
    output:
    tuple val(sample_name), path("${sample_name}.vcf.gz")
    
    script:
    """
      bgzip -c "${sample_name}".vcf > "${sample_name}".vcf.gz  

    """
}

process csiVCF_singularity {

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/bcftools:1.9--ha228f0b_4' :
    'quay.io/biocontainers/bcftools:1.9--ha228f0b_4' }"
    
    publishDir "${baseDir}/output/", mode: 'copy'

    input:
    tuple val(sample_name), path("${sample_name}.vcf.gz")

    output:
    tuple val(sample_name), path("${sample_name}.vcf.gz.csi")

    script:
    """
      bcftools index "${sample_name}".vcf.gz -o "${sample_name}".vcf.gz.csi    
    """
 }

process makeBcfConsensus_singularity {
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/bcftools:1.9--ha228f0b_4' :
    'quay.io/biocontainers/bcftools:1.9--ha228f0b_4' }"

    publishDir "${baseDir}/output/", mode: 'copy'
    
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

process makeGenomeCov_singularity {

   container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
   'https://depot.galaxyproject.org/singularity/bedtools:2.27.1--he941832_2' :
   'quay.io/biocontainers/bedtools:2.27.1--he941832_2' }"
      
   publishDir "${baseDir}/output/", mode: 'copy'
   
   input:
   tuple val(sample_name), path("${sample_name}.sorted.bam")
   path reference
   
   output:
   tuple val(sample_name), path("${sample_name}.bedGraph")
   
   script:
   """
   genomeCoverageBed -ibam "${sample_name}".sorted.bam -d -g "${reference}"/"${reference_name}".sizes > "${sample_name}".bedGraph
   """
   
 }

process makeCoverageMask_singularity {
   
   publishDir "${baseDir}/output/", mode: 'copy'
   
   input:
   tuple val(sample_name), path("${sample_name}.bedGraph")
   //path reference
   
   output:
   tuple val(sample_name), path("${sample_name}.bedGraph.FiveX")

   script:
   """
   awk '{ if (\$3 < 5) {print \$1"\t" \$2"\t" \$3"\t" "N" } else {print \$1"\t" \$2"\t" \$3"\t"}}' "${sample_name}".bedGraph > "${sample_name}".bedGraph.FiveX
   """
  
}

process maskWithNs_singularity {
   
   container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
   'https://depot.galaxyproject.org/singularity/seqtk:1.3--hed695b0_2' :
   'quay.io/biocontainers/seqtk:1.3--hed695b0_2' }"
   
   publishDir "${baseDir}/output/", mode: 'copy'
   
   input:
   tuple val(sample_name), path("${sample_name}.bedGraph.FiveX")
   tuple val(sample_name), path("${sample_name}.fasta")
   
   output:
   tuple val(sample_name), path("${sample_name}.fiveX.fasta")
      
   script:
   """
   seqtk mutfa "${sample_name}".fasta "${sample_name}".bedGraph.FiveX > "${sample_name}".fiveX.fasta 
   """
}

process queryVCF_singularity {

   container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
   'https://depot.galaxyproject.org/singularity/bcftools:1.9--h47928c2_1' :
   'quay.io/biocontainers/bcftools:1.9--h47928c2_1' }"
   
   publishDir "${baseDir}/output/", mode: 'copy'
   
   input:
   tuple val(sample_name), path("${sample_name}.vcf.gz")
   
   output:
   stdout
//   tuple val(sample_name), path("${sample_name}.snp.tsv")
   
   script:
   """
   bcftools query -f '''%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t[%AD]\t[%DP]\n''' "${sample_name}".vcf.gz 
   """
 //   | perl $PWD/Perl_scripts/bcftoolsQuery.pl > "${sample_name}".snp.tsv
}


workflow {
    
    Channel
        .fromFilePairs(params.inputPair, checkIfExists: true)
        .set { read_pairs_ch }
      
   // LOOKSY(read_pairs_ch, params.reference) | view
    
    mapResults = bowtie2map_singularity(read_pairs_ch, reference_path) 
    
    mapResults.view { "Bowtie2 Results: ${it}" }
    
    bamResults = sam2bam_singularity(mapResults)
    
    sortedBam = sortBam_singularity(bamResults)
    
    indexedBam = indexBam_singularity(sortedBam)
    
    myVCF = makeVCF_singularity(sortedBam, reference_path)   
    
    myVCF.view { "SNP calls: ${it}" }
    
    myVCFz = zipVCF_singularity(myVCF)

    myVCFcsi = csiVCF_singularity(myVCFz)

    fasta = makeBcfConsensus_singularity(myVCFz, myVCFcsi, reference_path)

    fasta.view { "Consensus FASTA: ${it}" }

    bedGr = makeGenomeCov_singularity(sortedBam, reference_path)

    bed5X = makeCoverageMask_singularity(bedGr)

    fastaN = maskWithNs_singularity(bed5X, fasta)

    //SNPs = queryVCF_singularity(myVCFz)

    queryVCF_singularity(myVCFz).view()
}
