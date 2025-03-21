#!/usr/bin/python

## Calls Singularity to install local bowtie2 as user-defined input or default bowtie2:2.5.1
## Also installs local perl:5.32, local samtools:1.9, and local pbgzip:2016.08.04

import sys, os, os.path, argparse, re, string, logging, warnings, time
import pathlib
import urllib.request
import subprocess

def check_website(url):  ## input is a website - website validation function by MS Copilot
    try:
        response = urllib.request.urlopen(url)
        if response.status == 200:
            print(f"The website, {url} exists.")
    except urllib.error.URLError as e:
        print(f"The website, {url} does not exist. Reason: {e.reason}")


def pull_singularity_image(image_url, output_path):  ## Spawns new process to pull singularity image
    try:
        # Construct the singularity pull command
        command = ['/usr/bin/singularity', 'pull', '--name', output_path, image_url]
        
        # Execute the command
        result = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        # Print the output and error (if any)
        print(result.stdout.decode())
        if result.stderr:
            print(result.stderr.decode())
        
        print(f"Image pulled successfully to {output_path}")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred: {e.stderr.decode()}")


logger = logging.getLogger("get_bowtie2_and_perl.py")
logger.setLevel(logging.INFO)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
#add formatter to ch
ch.setFormatter(formatter)
#add ch to logger
logger.addHandler(ch)

currentDir = os.getcwd()

parser = argparse.ArgumentParser(description='Pulls singularity container of bowtie2 v2.5 or later and builds a .sif called by bowtieConsensusFromPairedFastq.nf, learn_DSL2/map_bowtie_or_bwa.nf, and Python_scripts/prepReferenceAndReads.py', usage="python Singularity/get_bowtie2_and_perl.py --input 'https://depot.galaxyproject.org/singularity/bowtie2:2.5.4--he20e202_3' (default = no input, pull from https://depot.galaxyproject.org/singularity/bowtie2:2.5.1--py39h6fed5c7_1)")

parser.add_argument('--input', default="https://depot.galaxyproject.org/singularity/bowtie2:2.5.1--py39h6fed5c7_1")

args = parser.parse_args()

check_website(args.input)

logger.info("Begin pulling from Singularity repository, {}".format(args.input))

time.sleep(1)

#os.system("singularity pull {}".format(args.input))
bowtie_sif = 'my_bowtie2.sif'

try:
    with open(bowtie_sif, 'r+'):
        print("{} already built.".format(bowtie_sif))
except:
    pull_singularity_image("https://depot.galaxyproject.org/singularity/bowtie2:2.5.1--py39h6fed5c7_1", bowtie_sif)

time.sleep(1)        

#findContainer = args.input.split('/')
#cont = len(findContainer)
#logger.info("Built image from {}".format(findContainer[cont - 1]))
#os.system("singularity build my_bowtie2.sif {}".format(findContainer[cont - 1]))

logger.info("Successful build of bowtie2.sif")

time.sleep(1)

logger.info("Begin Singularity build of perl v5.32 and samtools v1.9")

check_website("https://depot.galaxyproject.org/singularity/perl:5.32")

perl_sif = 'my_perl5.32.sif'

try:
    with open(perl_sif, 'r+'):
        print("{} already built.".format(perl_sif))
except:
    pull_singularity_image("https://depot.galaxyproject.org/singularity/perl:5.32", perl_sif)

#os.system("singularity pull https://depot.galaxyproject.org/singularity/perl:5.32")

logger.info("Built image from https://depot.galaxyproject.org/singularity/perl:5.32")

#os.system("singularity build my_perl5.32.sif perl:5.32")

time.sleep(1)

check_website("https://depot.galaxyproject.org/singularity/samtools:1.9--h91753b0_8")

samtools_sif = 'samtools1.9.sif';

try:
    with open(samtools_sif, 'r+'):
        print("{} already built.".format(samtools_sif))
except:
    pull_singularity_image("https://depot.galaxyproject.org/singularity/samtools:1.9--h91753b0_8", samtools_sif)


#os.system("singularity pull https://depot.galaxyproject.org/singularity/samtools:1.9--h91753b0_8")
logger.info("Built image from https://depot.galaxyproject.org/singularity/samtools:1.9--h91753b0_8")
#os.system("singularity build my_samtools.sif samtools:1.9--h91753b0_8")

time.sleep(1)

check_website("https://depot.galaxyproject.org/singularity/htslib:1.19.1--h81da01d_2")

pbgzip_sif = 'htslib1.19.1.sif';

try:
    with open(pbgzip_sif, 'r+'):
        print("{} already built.".format(samtools_sif))
except:
    pull_singularity_image("https://depot.galaxyproject.org/singularity/htslib:1.19.1--h81da01d_2", pbgzip_sif)

time.sleep(1)

logger.info("Built image from https://depot.galaxyproject.org/singularity/htslib:1.19.1--h81da01d_2")



