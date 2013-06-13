#
# STUB not usable
#

#!/usr/bin/env python

import sys
import tempfile
import os

docstring= """
DESCRIPTION
    Split a file in multiple files on the bases of a column of
    identifiers (reverse of concatenate_bed.py)
    
    Input MUST be sorted by the identifying column. Use - to read from stdin

USAGE:
    deconcatenate.py <infile> <col-idx>

"""
if len(sys.argv) != 3:
    print(docstring)
    sys.exit(1)

if sys.argv[1] in ('-h', '--help'):
    print(docstring)
    sys.exit(0)

tmp= '/lustre/sblab/berald01/Tritume'## tempfile.mkdtemp(prefix= 'deconcatenate_')
outfiles= {}
if sys.argv[1] == '-':
    fin= sys.stdin
else:
    fin= open(sys.argv[1])

cfile= ''
start= True
while True:
    line= fin.readline()
    line= line.strip().split('\t')
    print(int(sys.argv[2]) - 1)
    id= line[int(sys.argv[2]) - 1]
    if id != cfile:
        if start:
            start= False
        else:
            fout.close()
        fname= os.path.join(tmp, id + '.txt')
        fout= open(fname, 'w')
        fout.write('\t'.join(line) + '\n')
        cfile= id
    else:
        fout.write('\t'.join(line) + '\n')
