#!/usr/bin/env python

docstring= """DESCRIPTION
    Prepare PE bam file for reading by R function bam2gr (package PING).
    Second mate reads get '/3' string appended to read name to comply with bam2gr.
    It is assumed that paired reads have exactly the same name. If reads have /1 and /2
    appended (output from bowtie/bismark) use cleanBamReadNames.py to remove
    /1 /2.

USAGE
    samtools view -h <in.bam> \
        | bam2grPrepareBam.py \
        | samtools view -S -b > <out.bam>

samtools view -h /lustre/sblab/berald01/repository/bam_clean/ear012_MNAse_10.bwa.mm9.test.bam \
    | bam2grPrepareBam.py \
    | samtools view -S -b - > ear012_MNAse_10.bwa.mm9.test2.bam
"""

import sys

fin= sys.stdin
for line in fin:
    line= line.rstrip('\n')
    if line.startswith('@'):
        print(line)
        continue
    line= line.split('\t')
    if int(line[1]) & 128 == 128:
        line[0]= line[0] + '/3'
    print('\t'.join(line))
        
sys.exit()