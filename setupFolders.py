#!/usr/bin/python

import sys, os, os.path, argparse, re, string, logging, warnings

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

parser = argparse.ArgumentParser(description='Creates the input folders, fileName/ and referenceFileName/, along with output folders, bowtieConsensInterm/ and bowtieConsensOutput/, expected as parameters for bowtieConsensusFromPairedFastq.nf', usage="python setupFolders.py --input filepath/fileNames.tar(.gz)")

parser.add_argument('--input', default='bowtieConsensTestFiles.tar', type=ext_check('.tar', 'tar.gz', argparse.FileType('r')))

## parse and print input files from .tar
args = parser.parse_args()
print(args.input.name)
extracted = re.sub(r'\.tar', '', args.input.name)

## declare input and output folders expected by bowtieConsensusFromPairedFastq.nf
fastqFolder = 'Noro_paired_files'
referenceFolder = 'reference_Noro_paired_files'
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

#os.system("cd {}".format(fastqFolder))
#os.system("gunzip {}/*{}".format(extracted, "fastq.gz"))
#os.system("cd {}".format(currentDir))

os.system("cp -v {}/*{} {}".format(extracted, "fasta", referenceFolder))

fastaFile = ''
files = os.listdir(referenceFolder)
for f in files:
        if(f.endswith('fasta')):
                fastaFile = f


response = os.popen("which bowtie2").read()
logger.info("Found bowtie2 in {}".format(response))
response = os.popen("which samtools").read()
logger.info("Found samtools in {}".format(response))
response = os.popen("which genomeCoverageBed").read()
logger.info("Found genomeCoverageBed in {}".format(response))
response = os.popen("which freebayes").read()
logger.info("Found freebayes in {}".format(response))

refBaseName = re.sub(r'\.fasta', '', fastaFile)

os.system("bowtie2-build {}/{} {}/{}".format(referenceFolder, fastaFile, referenceFolder, refBaseName))
logger.info("bowtie2 reference indexing complete.")

faidxName = fastaFile + ".fai"
sizesName = refBaseName + ".sizes"

os.system("samtools faidx {}/{}".format(referenceFolder, fastaFile))
os.system("cut -f 1,2 {}/{} > {}/{}".format(referenceFolder, faidxName, referenceFolder, sizesName))
logger.info("fasta reference length indexing complete.")

