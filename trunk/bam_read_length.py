#!/usr/bin/env python

import pysam
import sys

docstring="""DESCRIPTION
    Produce a histogram of read length from bam file
USAGE
    bam_read_length.py <sam|bam> <max-lines>    

    sam|bam: Input sam or bam file
    max-reads: Stop after reading this many reads having
        read length not None (default to read all the reads)

OUTPUT
    Tab separated file with: input file name, read length, read count.
"""

if len(sys.argv) not in (2, 3) or sys.argv[1] in ('-h', '--help'):
    sys.exit(docstring)

if len(sys.argv) == 3:
    max_reads= int(sys.argv[2])
else:
    max_reads= None

read_count= 0
read_len_hist= {}

samfile= pysam.Samfile(sys.argv[1], 'rb')
for line in samfile:
    rl= line.alen
    if rl is not None:
        read_count += 1
        read_len_hist[rl]= read_len_hist.get(rl, 0) + 1
        if max_reads is not None and read_count >= max_reads:
            break
    
lengths= sorted([k for k in read_len_hist])
for i in lengths:
    print('\t'.join([sys.argv[1], str(i), str(read_len_hist[i])]))

samfile.close()
sys.exit()