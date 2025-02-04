#!/usr/bin/python

import sys, os, os.path, argparse, re, string, logging, warnings, time
import pathlib
import urllib.request

def check_website(url):  ## input is a website - website validation function by MS Copilot
    try:
        response = urllib.request.urlopen(url)
        if response.status == 200:
            print(f"The website, {url} exists.")
    except urllib.error.URLError as e:
        print(f"The website, {url} does not exist. Reason: {e.reason}")


logger = logging.getLogger("get_bowtie2.py")
logger.setLevel(logging.INFO)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
#add formatter to ch
ch.setFormatter(formatter)
#add ch to logger
logger.addHandler(ch)

currentDir = os.getcwd()

parser = argparse.ArgumentParser(description='Pulls singularity container of bowtie2 v2.5 or later and builds a .sif called by bowtieConsensusFromPairedFastq.nf, learn_DSL2/map_bowtie_or_bwa.nf, and Python_scripts/prepReferenceAndReads.py', usage="python Singularity/get_bowtie2.py --input 'https://depot.galaxyproject.org/singularity/bowtie2:2.5.4--he20e202_3' (default = no input, pull from https://depot.galaxyproject.org/singularity/bowtie2:2.5.1--py39h6fed5c7_1)")

parser.add_argument('--input', default="https://depot.galaxyproject.org/singularity/bowtie2:2.5.1--py39h6fed5c7_1")

args = parser.parse_args()

check_website(args.input)

logger.info("Begin pulling from Singularity repository, {}".format(args.input))

time.sleep(1)

os.system("singularity pull {}".format(args.input))

time.sleep(1)

findContainer = args.input.split('/')

cont = len(findContainer)

logger.info("Build image from {}".format(findContainer[cont - 1]))

os.system("singularity build my_bowtie2.sif {}".format(findContainer[cont - 1]))

time.sleep(1)

logger.info("Successful build of bowtie2.sif")
