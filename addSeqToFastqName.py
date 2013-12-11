#!/usr/bin/env python

import sys
import gzip

docstring= """Add to the read name of a fastq file the first n and last m bases of each read
as barcode.

USAGE:
    addSeqToFastqName.py <in-fastq> <N-left> <M-right>
    
    in-fastq:
        Input fastq file or - to read from stdin
    N-left:
        Add this many bases from the left of the read to the read name
    M-right
        Add this many bases from the right of the read to the read name

OUTPUT FORMAT:
    The left and right barcodes go *after* the first blank space, separted by colon. E.g.
    
    @M00621:7:000000000-A6H18:1:1101:16086:1364:ACTG:ACTG 1:N:0:7

MEMO: Python code to extract the barcodes:
    rname= '@M00621:7:000000000-A6H18:1:1101:16086:1364:TATGATTATGGT:AGAGATCGGAAG 1:N:0:7'
    rname.split()[0].split(':')[-2:]
    >>> ['TATGATTATGGT', 'AGAGATCGGAAG']
"""

def addBarcodesToName(rname, seq, leftlen, rightlen):
    """Add barcode sequences to the read name *after* the first blank space.
    rname:
        Read name to edit (first line of fastq)
    seq:
        Sequence from which to extract barcode (2nd line fastq)
    leftlen, rightlen:
        Left and right length of the barcode.
    """
    if leftlen > len(seq):
        leftlen= len(seq)
    if rightlen > len(seq):
        rightlen= len(seq)

    xname= rname.split()
    prefix= xname[0]
    suffix= ' '.join(xname[1:])
    newname= prefix + ':' + seq[0:leftlen] + ':' + seq[(len(seq) - rightlen):]
    return(newname + ' ' + suffix)

# ------------------------------------------------------------------------------

if len(sys.argv) != 4 or sys.argv[1] in ('-h', '--help'):
    sys.exit(docstring)

try:
    leftlen= int(sys.argv[2])
    rightlen= int(sys.argv[3])
except ValueError:
    sys.exit('Integer expected as length of sequence to extract')

if int(sys.argv[2]) < 0 or int(sys.argv[3]) < 0:
    sys.exit('Need int >= 0')

if sys.argv[1] == '-':
    fin= sys.stdin
elif sys.argv[1].endswith('.gz'):
    fin= gzip.open(sys.argv[1])
else:
    fin= open(sys.argv[1])

i= 0
readList= []
for line in fin:
    line= line.strip()
    if i < 4:
        readList.append(line)
        i += 1
    else:
        readList[0]= addBarcodesToName(readList[0], readList[1], leftlen, rightlen)
        print('\n'.join(readList))   
        readList= [line]
        i= 1
sys.exit()