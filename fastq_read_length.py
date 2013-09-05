#!/usr/bin/env python

import sys
import gzip

docstring="""
DESCRIPTION
    Produce histrogram of read length from input fastq file (gzipped ok)

USAGE
    fastq_read_length.py <fastq> <max-reads>
    
    fastq: Input fastq file
    max-reads: Stop after reading this many reads (default to read all the reads)

OUTPUT
    Tab separated file with: input file name, read length, read count.
"""

if len(sys.argv) not in (2, 3) or sys.argv[1] in ('-h', '--help'):
    sys.exit(docstring)

if sys.argv[1].endswith('.gz'):
    fin= gzip.open(sys.argv[1])
else:
    fin= open(sys.argv[1])

if len(sys.argv) == 3:
    max_reads= int(sys.argv[2])
else:
    max_reads= None

n= -1
read_count= 0
read_len_hist= {}

for line in fin:
    if n % 4 == 0:
        read_count += 1
        rl= len(line.rstrip())
        read_len_hist[rl]= read_len_hist.get(rl, 0) + 1
        if max_reads is not None and read_count >= max_reads:
            break
    n += 1
    
lengths= sorted([k for k in read_len_hist])
for i in lengths:
    print('\t'.join([sys.argv[1], str(i), str(read_len_hist[i])]))

fin.close()
sys.exit()