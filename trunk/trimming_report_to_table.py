#!/home/berald01/.local/bin/python

import argparse
import os
import sys

parser = argparse.ArgumentParser(description= """
DESCRIPTION
    Convert the trimming report produced by cutadapt/trim_galore to a single
    line suitable for fitting in a database table.
    
EXAMPLE
    
DEPENDS-ON:

DEPENDS-ON-ME:

TODO:

""", formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('infile',
                   help='''Trimming report to parse.
                   ''')

args = parser.parse_args()

report= open(args.infile).readlines()
report= [x.strip('\n') for x in report]

print(report)

"""
SUMMARISING RUN PARAMETERS
==========================
Input filename: mjb001_bs_lmw.slx-5211.s_1_1.fq.gz
Quality Phred score cutoff: 20
Quality encoding type selected: ASCII+33
Adapter sequence: 'AGATCGGAAGAGC'
Minimum required adapter overlap (stringency): 1 bp
Minimum required sequence length for both reads before a sequence pair gets removed: 20 bp
Running FastQC on the data once trimming has completed
Running FastQC with the following extra arguments: --noextract --outdir /lustre/sblab/berald01/fastqc
Output file will be GZIP compressed


cutadapt version 1.0
Command line parameters: -f fastq -q 20 -O 1 -a AGATCGGAAGAGC mjb001_bs_lmw.slx-5211.s_1_1.fq.gz
Maximum error rate: 10.00%
   Processed reads: 38157
     Trimmed reads: 16305 ( 42.7%)
   Too short reads: 0 (  0.0% of processed reads)
    Too long reads: 0 (  0.0% of processed reads)
        Total time:      2.51 s
     Time per read:      0.07 ms

=== Adapter 1 ===

Adapter 'AGATCGGAAGAGC', length 13, was trimmed 16305 times.

Histogram of adapter lengths
length	count
1	7244
2	2165
3	591
4	206
5	34
6	27
7	19
8	21
9	55
10	14
11	21
12	3
13	5905


RUN STATISTICS FOR INPUT FILE: mjb001_bs_lmw.slx-5211.s_1_1.fq.gz
=============================================
38157 sequences processed in total
Sequences removed because they became shorter than the length cutoff of 20 bp:	0 (0.0%)
"""

sys.exit()
