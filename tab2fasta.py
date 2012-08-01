#!/usr/local/bin/python

docstring= """
DESCRIPTION
    Convert tabular to FASTA

USAGE:
    python tab2fasta.py <tab-file> <sequence column> <header column 1> <header column 2> <header column n>  > <outfile>
"""

import sys
if len(sys.argv) < 4:
    sys.exit('\nThree or more arguments required%s' %(docstring))
    
infile= open(sys.argv[1])
seqix= int(sys.argv[2]) - 1 
headerix= sys.argv[3:]
headerix= [(int(x) - 1) for x in headerix]

for line in infile:
    line= line.strip().split('\t')
    header= '>' + '_'.join([line[i] for i in headerix])
    print(header)
    print(line[seqix])

infile.close()
