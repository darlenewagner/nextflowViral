#!/usr/bin/python

import sys, os, os.path, argparse, re, string, logging, warnings, time
import pathlib

## Function: Checks existence of --refDir and --readDir
def readable_dir(prospective_dir):
        if not os.path.isdir(prospective_dir):
                raise argparse.ArgumentTypeError("readable_dir:{0} is not a valid path".format(prospective_dir))
        if os.access(prospective_dir, os.R_OK):
                if( not prospective_dir.endswith("/") ):
                        prospective_dir = prospective_dir + "/"
                return prospective_dir
        else:
                raise argparse.ArgumentTypeError("readable_dir:{0} is not a readable dir".format(prospective_dir))


## config file as optional input, default is 

def ext_check(expected_ext, openner):
        def extension(filename):
                if not (filename.lower().endswith(expected_ext)):
                        raise ValueError()
                return openner(filename)
        return extension

logger = logging.getLogger("prepReferenceOnly.py")
logger.setLevel(logging.INFO)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
#add formatter to ch
ch.setFormatter(formatter)
#add ch to logger
logger.addHandler(ch)

parser = argparse.ArgumentParser(description="Wrapper for reference indexing using samtools", usage="python prepReferenceOnly.py refDir/ ")

parser.add_argument('refDir', type=readable_dir)

args = parser.parse_args()

referenceFolder = args.refDir


logger.info("Reading folder of fasta reference {}".format(referenceFolder))

fasta = os.listdir(referenceFolder)

## Make sure there is only one fasta file in referenceFolder
confirmSingle = 0
fileToIndex = ''

for fas in fasta:
        if(re.search(r'(\.fasta|\.fas)?', fas)):
                print(fas)
                confirmSingle = confirmSingle + 1
        if(confirmSingle == 1):
                fileToIndex = fas

indexedFile = re.sub(r'\.fasta', '', fileToIndex)


if(confirmSingle > 1):
        logger.warning("Reference folder {} contains {} .fasta/.fas files. Only indexing {}.".format(referenceFolder, confirmSingle, indexedFile))

os.system("bowtie2-build {}/{} {}/{}".format(referenceFolder, fileToIndex, referenceFolder, indexedFile))
logger.info("bowtie2 reference indexing complete.")

faidxName = fileToIndex + ".fai"
sizesName = indexedFile + ".sizes"

os.system("samtools faidx {}/{}".format(referenceFolder, fileToIndex))
os.system("cut -f 1,2 {}/{} > {}/{}".format(referenceFolder, faidxName, referenceFolder, sizesName))
logger.info("fasta reference length indexing complete.")

