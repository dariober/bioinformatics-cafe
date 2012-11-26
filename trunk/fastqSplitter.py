#!/usr/bin/env python

import sys
import gzip

docstring= """\nDESCRIPTION
    Split the sequences and qualty strings of a fastq files in chunks of no less
    then a defined length. Input file can be gzipped.
    NB: Comments are stripped
        
USAGE
    fastqSplitter.py <size> <fastqfile> > <out fastq>

ARGUMENTS:  
    size:
        Each read will be split in chunks of this size or larger.
    fastqfile:
        Input fastq. Can be read as gzip.
        
EXAMPLE
    echo "@seq-1" > test.fq
    echo "AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTTTTTTT" >> test.fq
    echo "+comment" >> test.fq
    echo "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII" >> test.fq
    echo "@seq-2" >> test.fq
    echo "NAAAAAAAAANCCCCCCCCCNGGGGGGGGGTTTTTTTTTT" >> test.fq
    echo "+comment" >> test.fq
    echo "########################################" >> test.fq
    echo "@seq-3" >> test.fq
    echo "NAAAAAA" >> test.fq
    echo "+comment" >> test.fq
    echo "#######" >> test.fq

fastqSplitter.py 10 test.fq > test.split.fq
fastqSplitter.py 10 test.fq | gzip > test.split.fq.gz

"""

if sys.argv[1] in ('-h', '--help'):
    sys.exit(docstring)

if len(sys.argv) != 3:
    print('\n*** Error ***\n')
    sys.exit(docstring)

# -----------------------------------------------------------------------------

class Fastq:
    def __int__(self):
        name= None
        sequence= None
        comment= None
        quality= None
    def fqPrint(self):
        print(self.name)
        print(self.sequence)
        print(self.comment)
        print(self.quality)

def readFastqLine(fqin):
    """Read 4 lines from fastq file handle fqin and returns it as Fastq object or
    None if EOF is reached"""
    fq= Fastq()
    fq.name= fqin.readline()
    if fq.name == '':
        return(None)
    else:
        fq.name= fq.name.strip()
    fq.sequence= fqin.readline().strip()
    fq.comment= fqin.readline().strip()
    fq.quality= fqin.readline().strip()
    return(fq)

    
def getSplittingPoints(x, splitlen):
    """Split string x in chunks of length 'splitlen' (int). Last chunk is rounded
    in excess so it will be of length: splitlen <= last < 2*splitlen.
    Returns:
        List of ints where splitting will occur
    """
    l= len(x)
    splits= [s * splitlen for s in range(1, l/splitlen)]
    splits.append(l)
    return(splits)
    
    
def splitter(fq, splitlen):
    """Split sequences and qualities in fq (Fastq obj) in chunks of size splitlen
    (int) or larger
    Returns:
        List of Fastq objects, one for each splitted sequnces
    """
    fqs= []
    splittingpoints= getSplittingPoints(fq.sequence, splitlen)
    startpoints= [0] + splittingpoints[0:-1]
    for s, e in zip(startpoints, splittingpoints):
        fastq= Fastq()
        fastq.name= fq.name + '-' + str(s)
        fastq.sequence= fq.sequence[s:e]
        fastq.comment= '+'
        fastq.quality= fq.quality[s:e]
        fqs.append(fastq)
    return(fqs)
# ------------------------------------------------------------------------------

splitlen= int(sys.argv[1])
fastq= sys.argv[2]

if fastq.endswith('.gz'):
    fqin= gzip.open(fastq)
else:
    fqin= open(fastq)

while True:
    fq= readFastqLine(fqin)
    if fq is None:
        break
    splitFastqList= splitter(fq, splitlen)
    for q in splitFastqList:
        q.fqPrint()

fqin.close()
