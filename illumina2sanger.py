#!/usr/local/bin/python

docstring= """
SYNOPSIS

    illumina2sanger.py <illumina-fastq> > <sanger-fastq>

DESCRIPTION

    Convert Illumina >1.3 and <1.9 encoding to Sanger.
    Throw an exception if characters less than ascii 64 (@) or above 126 (~) are found.
    
    If input file is gzip'd it will be decompressed and recompressed.

EXAMPLES

    illumina2sanger.py illumina.fq > sanger.fq
    illumina2sanger.py illumina.fq.gz > sanger.fq           <<< Input zipped, output unzipped
    illumina2sanger.py illumina.fq.gz | gzip > sanger.fq.gz <<< Input and output zipped

MEMO from Wikipedia:

    Encoding

    * Sanger format can encode a Phred quality score from 0 to 93 using ASCII 33 to 126
    * Solexa/Illumina 1.0 format can encode a Solexa/Illumina quality score from -5 to 62 using ASCII 59 to 126
    * Illumina 1.3+ format can encode a Phred quality score from 0 to 62 using ASCII 64 to 126.
    * The Phred scores 0 to 2 in Illumina 1.5+ have a slightly different meaning. 

"""

import sys
import subprocess
import os

if len(sys.argv) != 2 or sys.argv[1] == '-h' or sys.argv[1] == '--help':
    print(docstring)
    sys.exit()
    
def illumina2sanger(base):
    return(chr(ord(base) - 31))

infastq_file= sys.argv[1]
if infastq_file.endswith('.gz'):
    p= subprocess.Popen('gunzip %s' %(sys.argv[1]), shell= True)
    p.wait()
    fastq_file= infastq_file[:-3]
else:
    fastq_file= sys.argv[1]
    
fastq= open(fastq_file)

while True:
    qname= fastq.readline().rstrip('\n\r')
    if qname == '':
        break
    print(qname)
    print(fastq.readline().rstrip('\n\r'))
    print(fastq.readline().rstrip('\n\r'))
    qual= fastq.readline().rstrip('\n\r')
    qual_sanger= [illumina2sanger(x) for x in qual if ord(x) >= 64 and ord(x) <= 126]
    if len(qual_sanger) != len(qual):
        sys.exit('Invalid Illumina encoding for line %s' %(''.join(qual)))
    qual_sanger= ''.join(qual_sanger)
    print(qual_sanger)
fastq.close()

if sys.argv[1].endswith('.gz'):
    p= subprocess.Popen('gzip %s' %(fastq_file), shell= True)
    p.wait()

sys.exit()