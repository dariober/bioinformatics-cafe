#!/usr/bin/env python

import sys
import pysam

docstring= """Filter sam file (stream) to get the reads name listed in txt
    USAGE
        getReadsByName.py <bam> <reads.txt>

    reads.txt is file where first column is the read name, other columns are ignored
    Output is bam file printed to stdout.
    
    TO TEST:
        samtools view mybam.bam | head -n 10 > test.sam
        getReadsByName.py <mybam.bam> <test.sam>
"""

if len(sys.argv) != 3:
    sys.exit(docstring)

if sys.argv[1].endswith('.sam'):
    samfile = pysam.Samfile(sys.argv[1], "r")
elif sys.argv[1].endswith('.bam'):
    samfile = pysam.Samfile(sys.argv[1], "rb")
elif sys.argv[1] == '-':
    samfile = pysam.Samfile("-", "r")
else:
    sys.exit('Extension of input file must be .sam or .bam')

outfile = pysam.Samfile( "-", "wb", template = samfile)

#if args.outbam:
#    outfile = pysam.Samfile( "-", "wb", template = samfile)
#else:
#    outfile = pysam.Samfile( "-", "w", template = samfile)

## Read in file of read names
readNames= []
fin= open(sys.argv[2])
for line in fin:
    line= line.strip().split('\t')[0]
    readNames.append(line)
fin.close()
readNames= tuple(readNames)

for read in samfile:
    name= read.qname
    if name in readNames:
        outfile.write(read)
samfile.close()
outfile.close()
sys.exit()
