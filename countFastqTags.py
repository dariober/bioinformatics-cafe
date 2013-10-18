#!/usr/bin/env python

docstring= """
DESCRIPTION
    Count the occurences of each sequence in a fastq file.

USAGE
    countFastqTags.py <fastq>

OUTPUT
    To stdout a table of sequence and count tab separated
"""

import argparse
import sys
import os
import gzip

if len(sys.argv) != 2 or sys.argv[1] in ('-h', '--help', '-help'):
    sys.exit(docstring)

fastq= sys.argv[1]

if fastq.endswith('.gz'):
    fin= gzip.open(fastq)
else:
    fin= open(fastq)

seqdict= {}
i= -1
for line in fin:
    if i % 4 == 0:
        line= line.strip()
        if line in seqdict:
            seqdict[line] += 1
        else:
            seqdict[line]= 1
    i += 1

fin.close()

for k in sorted(seqdict.keys()):
    print(k + '\t' + str(seqdict[k]))

sys.exit()