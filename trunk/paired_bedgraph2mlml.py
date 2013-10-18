#!/usr/bin/env python

import sys

docstring= """DESCRIPTION
    Cut & Reformat bed file to produce two input files for mlml.

INPUT
    File obtained by pasting side by side two bedgraph files from mpileup2methylation.py.
    Example:
    chr1    14653   14654   100.0   1       1       -       chr1    14653   14654   0.0     0       1       -
    chr1    14699   14700   100.0   2       2       -       chr1    14699   14700   100.0   1       1       -
USAGE
    paired_bedgraph2mlml.py <infile> <outfile-1> <outfile-2>
"""

if len(sys.argv) != 4:
    sys.exit(docstring)
    
fin= open(sys.argv[1])
fout_1= open(sys.argv[2], 'w')
fout_2= open(sys.argv[3], 'w')

for line in fin:
    line= line.strip().split('\t')
    outline_1= [line[0], line[1], line[2], 'CpG:' + line[5], str(float(line[4]) / float(line[5])), line[6]]
    outline_2= [line[0], line[1], line[2], 'CpG:' + line[12], str(float(line[11]) / float(line[12])), line[13]]
    fout_1.write('\t'.join(outline_1) + '\n')
    fout_2.write('\t'.join(outline_2) + '\n')

fin.close()
fout_1.close()
fout_2.close()

sys.exit()