#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/* Usage
     Preferred: nextflow run consensHiCovFromPairedFastq.singularity.nf --reference <fasta file> --inputPair <paired fastq files>
     - For web-connection-dependent Singularity containerization which calls required container separately within each process definition
     
     Expert: nextflow run consensHiCovFromPairedFastq.singularity.nf --local --reference <fasta file> --inputPair <paired fastq files>
     - Use flag '--local' to run without pulling and building from Singularity, only recommended when prerequisite programs are available locally
 */

params.local = false 

params.reference = "${baseDir}/bowtieConsensTestFiles/eng_live_atten_poliovirus/MZ245455.1"
reference_name = file(params.reference).name
reference_idx = "${reference_name}.fai"
reference_path = file(params.reference).parent

params.inputPair = "${baseDir}/bowtieConsensTestFiles/eng_live_atten_poliovirus/polio-sample-8_S13_R{1,2}_001.fastq.gz"
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


// -- Singularity-Driven Processes --
// REQUIREMENTS: Only 3 software packages.
// Singularity version 1.4.0-1.el8 or higher
// nextflow v24.04.2 or higher
// Plus, your favorite version of perl

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
 //  stdout
   tuple val(sample_name), path("${sample_name}.messy.snp.tsv")
     
   
   script:
   """
   bcftools query -f '''%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t[%AD]\t[%DP]\n''' "${sample_name}".vcf.gz -o "${sample_name}".messy.snp.tsv
   """
 //   | perl $PWD/Perl_scripts/bcftoolsQuery.pl > "${sample_name}".snp.tsv
}


process callPerl {

   publishDir "${baseDir}/output/", mode: 'copy'

   input:
   tuple val(sample_name), path("${sample_name}.messy.snp.tsv")
   
   output:
   tuple val(sample_name), path("${sample_name}.snp.tsv")

   script:
   """
   cat "${sample_name}".messy.snp.tsv | perl $PWD/Perl_scripts/bcftoolsQuery.pl > "${sample_name}".snp.tsv
   """
}


// -- Locally-Driven Processes --
// WARNING: Local installation of the following prerequisites/dependencies will be validated during runtime:
// nextflow v24.04.2 or higher
// bowtie2/2.3.5.1 or higher
// samtools/1.9
// bcftools/1.9
// htslib/1.19.1 or higher
// bedtools/2.27.1
// seqtk/1.3
// plus, your favorite version of perl




process bowtie2map_local {

    publishDir "${baseDir}/intermediate/", mode: 'copy'
    
    input:
    tuple val(sample_name), path(reads)
    path reference
    val foundIt

    output:
    tuple val(sample_name), path("${sample_name}.sam")    

    when:
    foundIt.contains("true")

    script:
    """
    bowtie2 --no-unal --no-mixed -x "${reference}"/"${reference_name}" -1 "${reads[0]}" -2 "${reads[1]}" > "${sample_name}.sam"
    """
}

process sam2bam_local {

    publishDir "${baseDir}/intermediate/", mode: 'copy'
    
    input:
    tuple val(sample_name), path("${sample_name}.sam")
    val foundIt

    output:
    tuple val(sample_name), path("${sample_name}.bam")

    when:
    foundIt.contains("true")
    
    script:
    """
    samtools view "${sample_name}".sam -o "${sample_name}".bam
    """

}

process sortBam_local {

    publishDir "${baseDir}/intermediate/", mode: 'copy'

    input:
    tuple val(sample_name), path("${sample_name}.bam")
    val foundIt
    
    output:
    tuple val(sample_name), path("${sample_name}.sorted.bam")
    
    when:
    foundIt.contains("true")

    script:
    """
    samtools sort "${sample_name}".bam -o "${sample_name}".sorted.bam
    """
}

process indexBam_local {
    
    publishDir "${baseDir}/intermediate/", mode: 'copy'
    
    input:
    tuple val(sample_name), path("${sample_name}.sorted.bam")
    val foundIt
    
    output:
    tuple val(sample_name), path("${sample_name}.sorted.bam.bai")

    when:
    foundIt.contains("true")
    
    script:
    """
    samtools index -b "${sample_name}".sorted.bam "${sample_name}".sorted.bam.bai
    """

}


process makeVCF_local {

    publishDir "${baseDir}/output/", mode: 'copy'

    input:
    tuple val(sample_name), path("${sample_name}.sorted.bam")
    path reference
    val foundIt
    
    output:
    tuple val(sample_name), path("${sample_name}.vcf")

    when:
    foundIt.contains("true")

    script:
    """
    bcftools mpileup -d 35000 -Ob -f "${reference}"/"${reference_name}".fasta -Q 20 -q 20 --annotate FORMAT/AD,FORMAT/DP "${sample_name}".sorted.bam | bcftools call -mv --ploidy 1 -Ov --output "${sample_name}".vcf
    """
    
}

process zipVCF_local {
    
    publishDir "${baseDir}/output/", mode: 'copy'
    
    input:
    tuple val(sample_name), path("${sample_name}.vcf")
    val foundIt
    
    output:
    tuple val(sample_name), path("${sample_name}.vcf.gz")

    when:
    foundIt.contains("true")

    script:
    """
      bgzip -c "${sample_name}".vcf > "${sample_name}".vcf.gz  

    """
}

process csiVCF_local {

    publishDir "${baseDir}/output/", mode: 'copy'

    input:
    tuple val(sample_name), path("${sample_name}.vcf.gz")
    val foundIt

    output:
    tuple val(sample_name), path("${sample_name}.vcf.gz.csi")

    when:
    foundIt.contains("true")

    script:
    """
      bcftools index "${sample_name}".vcf.gz -o "${sample_name}".vcf.gz.csi    
    """
 }


process makeBcfConsensus_local {
    
    publishDir "${baseDir}/output/", mode: 'copy'
    
    input:
    tuple val(sample_name), path("${sample_name}.vcf.gz")
    tuple val(sample_name), path("${sample_name}.vcf.gz.csi")
    path reference
    val foundIt
    
    output:
    tuple val(sample_name), path("${sample_name}.fasta")

    when:
    foundIt.contains("true")

    script:
    """
    cat "${reference}"/"${reference_name}".fasta | bcftools consensus "${sample_name}".vcf.gz --sample "${sample_name}".sorted.bam -o "${sample_name}".fasta
    """
    
}

process makeGenomeCov_local {

    publishDir "${baseDir}/output/", mode: 'copy'
   
    input:
    tuple val(sample_name), path("${sample_name}.sorted.bam")
    path reference
    val foundIt
   
    output:
    tuple val(sample_name), path("${sample_name}.bedGraph")

    when:
    foundIt.contains("true")

    script:
    """
    genomeCoverageBed -ibam "${sample_name}".sorted.bam -d -g "${reference}"/"${reference_name}".sizes > "${sample_name}".bedGraph
    """
   
 }

process makeCoverageMask_local {
   
    publishDir "${baseDir}/output/", mode: 'copy'
    
    input:
    tuple val(sample_name), path("${sample_name}.bedGraph")
    val foundIt
    //path reference
   
    output:
    tuple val(sample_name), path("${sample_name}.bedGraph.FiveX")

    when:
    foundIt.contains("true")

    script:
    """
    awk '{ if (\$3 < 5) {print \$1"\t" \$2"\t" \$3"\t" "N" } else {print \$1"\t" \$2"\t" \$3"\t"}}' "${sample_name}".bedGraph > "${sample_name}".bedGraph.FiveX
    """
  
}

process maskWithNs_local {
   
    publishDir "${baseDir}/output/", mode: 'copy'
   
    input:
    tuple val(sample_name), path("${sample_name}.bedGraph.FiveX")
    tuple val(sample_name), path("${sample_name}.fasta")
    val foundIt
   
    output:
    tuple val(sample_name), path("${sample_name}.fiveX.fasta")

    when:
    foundIt.contains("true")

    script:
    """
    seqtk mutfa "${sample_name}".fasta "${sample_name}".bedGraph.FiveX > "${sample_name}".fiveX.fasta 
    """
}

process queryVCF_local {

    publishDir "${baseDir}/output/", mode: 'copy'
   
    input:
    tuple val(sample_name), path("${sample_name}.vcf.gz")
    val foundIt

    output:
    //  stdout
    tuple val(sample_name), path("${sample_name}.messy.snp.tsv")
     
    when:
    foundIt.contains("true")

    script:
    """
    bcftools query -f '''%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t[%AD]\t[%DP]\n''' "${sample_name}".vcf.gz -o "${sample_name}".messy.snp.tsv
    """
    //   | perl $PWD/Perl_scripts/bcftoolsQuery.pl > "${sample_name}".snp.tsv
}



include { checkExecutables as checkExecutables0 } from './modules/reusable'
include { checkExecutables as checkExecutables1 } from './modules/reusable'
include { checkExecutables as checkExecutables2 } from './modules/reusable'
include { checkExecutables as checkExecutables3 } from './modules/reusable'
include { checkExecutables as checkExecutables4 } from './modules/reusable'
include { checkExecutables as checkExecutables5 } from './modules/reusable'
include { checkExecutables as checkExecutables6 } from './modules/reusable'

workflow {
    
    Channel
        .fromFilePairs(params.inputPair, checkIfExists: true)
        .set { read_pairs_ch }
      
   // LOOKSY(read_pairs_ch, params.reference) | view
    
    if(params.local)
       {
         // validate bowtie2 installation
         foundIt0 = checkExecutables0( 'bowtie2' )
         foundIt0.view { "bowtie2 executable found: ${it}" }
    
         mapResults = bowtie2map_local(read_pairs_ch, reference_path, foundIt0) 
         mapResults.view { "Bowtie2 Results: ${it}" }

         // validate samtools installation
         foundIt1 = checkExecutables1( 'samtools' )
         foundIt1.view { "samtools executable found: ${it}" }
	 
         bamResults = sam2bam_local(mapResults, foundIt1)
	 sortedBam = sortBam_local(bamResults, foundIt1)
	 indexedBam = indexBam_local(sortedBam, foundIt1)

         // validate bcftools installation
         foundIt2 = checkExecutables2( 'bcftools' )
         foundIt2.view { "bcftools executable found: ${it}" }
	 
         myVCF = makeVCF_local(sortedBam, reference_path, foundIt2)   
         myVCF.view { "SNP calls: ${it}" }

         // validate htslib installation
         foundIt3 = checkExecutables3( 'bgzip' )
         foundIt3.view { "bgzip executable found: ${it}" }

         myVCFz = zipVCF_local(myVCF, foundIt3)
         myVCFcsi = csiVCF_local(myVCFz, foundIt3)

         // validate bedtools installation
         foundIt4 = checkExecutables4( 'genomeCoverageBed' )
         foundIt4.view { "genomeCoverageBed executable found: ${it}" }

         fasta = makeBcfConsensus_local(myVCFz, myVCFcsi, reference_path, foundIt4)
         fasta.view { "Consensus FASTA: ${it}" }
         bedGr = makeGenomeCov_local(sortedBam, reference_path, foundIt4)

         // validate seqtk installation
         foundIt5 = checkExecutables5( 'seqtk' )
         foundIt5.view { "seqtk executable found: ${it}" }

         bed5X = makeCoverageMask_local(bedGr, foundIt5)
         fastaN = maskWithNs_local(bed5X, fasta, foundIt5)

         // validate perl installation
         foundIt6 = checkExecutables6( 'perl' )
         foundIt6.view { "perl executable found: ${it}" }
         
         SNPs = queryVCF_local(myVCFz, foundIt6)

      }
      else
      {
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
         SNPs = queryVCF_singularity(myVCFz)
	
      }	

    filtered = callPerl(SNPs)

    //queryVCF_singularity(myVCFz).view()
}
