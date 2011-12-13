#!/usr/bin/python

""" Convert a list of sam files to BAM using samtools command
    "samtools view -b -S -o myoutput.bam myinput.sam"
    SAM files have header. (if not use "samtools faidx my_reference.fasta")

"""

## Path to samtools
samdir= '/exports/work/vet_roslin_nextgen/dario/samtools/samtools-0.1.10'

## List of sam files to convert
samfiles= ['LN_109.fq.sam', 'LN114.fq.sam', 'LN124.fq.sam', 'LN_130.fq.sam', 'LN_146.fq.sam', 'LN_173.fq.sam', 'LN_183.fq.sam', 'LN_20B.fq.sam', 'LN_21.fq.sam', 'LN_38.fq.sam', 'LN_47.fq.sam', 'LN_50.fq.sam', 'LN_57.fq.sam', 'LN_58.fq.sam', 'LN_92.fq.sam']

## Output is the same as input with .sam replaced by .bam

# -----------------------------------------------------------------------------

import os
import re

samtools= os.path.join(samdir, 'samtools view -b -S -o ')

for sam in samfiles:
    bam= re.sub('.sam$', '.bam', sam)
    cmd= samtools + bam + ' ' + sam
    print('Executing: ' + cmd)
    os.system(cmd)    

