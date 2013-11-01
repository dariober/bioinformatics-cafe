#!/usr/bin/env python

import sys
import zlib

docstring= """Read a file of sequences, one per line and return sequence
complexity
USAGE:
    sequenceComplexity.py <infile>
Use - to read from stdin.
OUTPUT
    <sequence> <sequence length> <compressed seq len> <seq len / compressed seq len>
SEE ALSO:
    http://www.biostars.org/p/44545/#44572
"""

def seqComp(s):
    sl= len(s)
    cl= len(zlib.compress(s))
    cr= cl / float(sl)
    return([s, sl, cl, cr])

if len(sys.argv) != 2 or sys.argv[1] in ['-h', '--help']:
    sys.exit(docstring)

if sys.argv[1] == '-':
    fin= sys.stdin
else:
    fin= open(sys.argv[1])

for line in fin:
    s= line.strip().split('\t')[0]
    seqcomp= seqComp(s)
    print('\t'.join([str(x) for x in seqcomp]))

sys.exit()