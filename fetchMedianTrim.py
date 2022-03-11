#! /apps/x86_64/python/3.9.1/bin/python

import sys, csv, statistics
from itertools import islice


head = []

## First input file is a blastn output -outfmt 6
with open(sys.argv[1], mode='r') as myBlast:
    head = list(islice(myBlast, 10))

fivePrim = []
threePrim = []

topFivePrim = 0
topThreePrim = 0

for line in head:
    findCoords = line.split('\t')
    if((topFivePrim == 0) and (topThreePrim ==0)):
        topFivePrim = int(findCoords[6])
        topThreePrim = int(findCoords[7])
    if((int(findCoords[7]) - int(findCoords[6]) + 10) >= (topThreePrim - topFivePrim)):
        fivePrim.append(int(findCoords[6]))
        threePrim.append(int(findCoords[7]))


startCoord = int(statistics.median(fivePrim))
endCoord = int(statistics.median(threePrim))

header = ''
sequence = []

## Second input is fasta-formatted nucleotide file
with open(sys.argv[2], mode='r') as myFasta:
    header = myFasta.readline()
    for line in myFasta:
        line.strip()
        sequence.append(line.strip())

print(header, end="")
fullSeq = ''.join(sequence)
trimSeq = ''

print(fullSeq[startCoord:endCoord])



