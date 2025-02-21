#!/usr/bin/python

import sys, os, os.path, argparse, re, string, logging, warnings, time
import pathlib

## one input parameter expected, default is 

def ext_check(expected_ext, expected_ext3, openner):
        def extension(filename):
                if not (filename.lower().endswith(expected_ext) or filename.lower().endswith(expected_ext3)):
                        raise ValueError()
                return openner(filename)
        return extension

logger = logging.getLogger("setupFolders.py")
logger.setLevel(logging.INFO)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
#add formatter to ch
ch.setFormatter(formatter)
#add ch to logger
logger.addHandler(ch)

currentDir = os.getcwd()

parser = argparse.ArgumentParser(description='Creates the input folders under bowtieConsensTestFiles/, along with output folders, intermediate/ and output/, expected as parameters for bowtieConsensusFromPairedFastq.nf', usage="python setupFolders.py --input filepath/fileNames.tar(.gz)")

parser.add_argument('--input', default='bowtieConsensTestFiles.tar', type=ext_check('.tar', 'tar.gz', argparse.FileType('r')))

## parse and print input files from .tar
args = parser.parse_args()
print(args.input.name)
extracted = re.sub(r'\.tar', '', args.input.name)

## declare input and output folders expected by bowtieConsensusFromPairedFastq.nf
fastqFolder = 'D70_paired_files'
referenceFolder = 'D70_reference_genome'
intermFolder = 'bowtieConsensInterm'
outputFolder = 'bowtieConsensOutput'

logger.info("Creating input folders {} and {}".format(fastqFolder, referenceFolder))
os.system("mkdir {}".format(fastqFolder))
os.system("mkdir {}".format(referenceFolder))

logger.info("Creating input folders {} and {}".format(intermFolder, outputFolder))
os.system("mkdir {}".format(intermFolder))
os.system("mkdir {}".format(outputFolder))

logger.info("Opening data archive {}".format(args.input.name))
os.system("tar -xvf {}".format(args.input.name))
os.system("cp -v {}/*{} {}".format(extracted, "fastq.gz", fastqFolder))

logger.info("Gunzipping .fastq.gz in {}".format(fastqFolder))
files = os.listdir(fastqFolder)
for f in files:
        if(f.endswith('gz')):
                newName = re.sub(r'\.gz', '', f)
                os.system('gunzip -c {}/{} > {}/{}'.format(fastqFolder, f, fastqFolder, newName))
                os.system('rm -v {}/{}'.format(fastqFolder, f))

os.system("cp -v {}/*{} {}".format(extracted, "fasta", referenceFolder))

fastaFile = ''
files = os.listdir(referenceFolder)
for f in files:
        if(f.endswith('fasta')):
                fastaFile = f

#if(re.search(r'scicomp\/home-pure', currentDir)):
#        logger.info("Attempting installation of bowtie2 from {}".format(currentDir))
#        os.system("module load bowtie2/2.3.5.1")

response = os.popen("which bowtie2").read()
myPath = os.path.split(response)

haltingVariable = False
packages = { 'bowtie2' : False, 'samtools' : False, 'genomeCoverageBed' : False, 'freebayes' : False, 'nextflow' : False }

## The next four blocks check for presence of prerequsite packages ##

## Check Bowtie2
if not os.path.exists(myPath[0]):
        logger.warning("Package bowtie2 not installed and/or not found in $PATH!!!")
        logger.warning("If home directory is {}, try module load bowtie2/2.3.5.1".format(currentDir))
        haltingVariable = True
else:
        logger.info("Found bowtie2 in {}".format(response))
        packages['bowtie2'] = True

response = os.popen("which samtools").read()
myPath = response.strip()

## Check samtools
if not os.path.isfile(myPath):
        logger.warning("Package samtools not installed and/or not found in $PATH!!!")
        haltingVariable = True
else:
        logger.info("Found samtools in {}".format(response))
        packages['samtools'] = True

## Check genomeCoverageBed
response = os.popen("which genomeCoverageBed").read()
myPath = response.strip()
if not os.path.isfile(myPath):
        logger.warning("Package genomeCoverageBed not installed and/or not found in $PATH!!!")
        haltingVariable = True
else:
        logger.info("Found genomeCoverageBed in {}".format(response))
        packages['genomeCoverageBed'] = True

## Check freebayes
response = os.popen("which freebayes").read()
myPath = response.strip()
if not os.path.isfile(myPath):
        logger.warning("Package freebayes not installed and/or not found in $PATH!!!")
        haltingVariable = True
else:
        logger.info("Found freebayes in {}".format(response))
        packages['freebayes'] = True

## Check nextflow
response = os.popen("which nextflow").read()
myPath = response.strip()
if not os.path.isfile(myPath):
        logger.warning("Package nextflow not installed and/or not found in $PATH!!!")
        haltingVariable = True
else:
        logger.info("Found nextflow in {}".format(response))
        packages['nextflow'] = True

## Terminate setup if one or more prerequistes are missing ##
if(haltingVariable == True):
        logger.error("The following packages are required to complete setup:")
        time.sleep(0.5)
        program = list(packages.keys())
        for p in program:
                if(packages[p] == False):
                        print("\t" + p)
        time.sleep(0.5)
        logger.error("Setup for bowtieConsensusFromPairedFastq.nf terminating unsuccessfully...")
        time.sleep(0.5)
        sys.exit(1)

refBaseName = re.sub(r'\.fasta', '', fastaFile)

os.system("bowtie2-build {}/{} {}/{}".format(referenceFolder, fastaFile, referenceFolder, refBaseName))
logger.info("bowtie2 reference indexing complete.")

faidxName = fastaFile + ".fai"
sizesName = refBaseName + ".sizes"

os.system("samtools faidx {}/{}".format(referenceFolder, fastaFile))
os.system("cut -f 1,2 {}/{} > {}/{}".format(referenceFolder, faidxName, referenceFolder, sizesName))
logger.info("fasta reference length indexing complete.")

