#!/usr/local/bin/python

docstring= """
SYNOPSIS

    illumina2sanger.py <illumina-fastq or stdin> > <sanger-fastq>

DESCRIPTION

    Convert Illumina >1.3 and <1.9 encoding to Sanger.
    Throw an exception if characters less than ascii 64 (@) or above 126 (~) are found.
    Use - to read input file from stream
    Note: Comment lines are rejected if identical to header.
    
EXAMPLES

    illumina2sanger.py illumina.fq > sanger.fq
    gunzip -c illumina.fq.gz | illumina2sanger.py - > sanger.fq           <<< Input zipped, output unzipped
    gunzip -c | illumina2sanger.py - | gzip > sanger.fq.gz <<< Input and output zipped

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
import gzip

if len(sys.argv) != 2 or sys.argv[1] == '-h' or sys.argv[1] == '--help':
    print(docstring)
    sys.exit()
    
def illumina2sanger(base):
    return(chr(ord(base) - 31))

fastq_file= sys.argv[1]

delete_unzipped= False
"""
if infastq_file.endswith('.gz'):
    fastq_file= infastq_file[:-3]
    p= subprocess.Popen('gunzip -c %s > %s' %(sys.argv[1], fastq_file), shell= True)
    p.wait()
    delete_unzipped= True
else:
    fastq_file= sys.argv[1]
"""

if fastq_file == '-':
    fastq= sys.stdin
elif fastq_file.endswith('.gz'):
    fastq= gzip.open(fastq_file)
else:
    fastq= open(fastq_file)

while True:
    qname= fastq.readline().rstrip('\n\r')
    if qname == '':
        break
    seqline= fastq.readline().rstrip('\n\r')
    cmt_line= fastq.readline().rstrip('\n\r')
    if cmt_line[1:] == qname[1:]:
        cmt_line= '+'
    qual= fastq.readline().rstrip('\n\r') 
    qual_sanger= [illumina2sanger(x) for x in qual if ord(x) >= 64 and ord(x) <= 126]
    if len(qual_sanger) != len(qual):
        sys.exit('Invalid Illumina encoding for line %s' %(''.join(qual)))
    qual_sanger= ''.join(qual_sanger)
    print(qname)
    print(seqline)
    print(cmt_line)
    print(qual_sanger)
fastq.close()

#if delete_unzipped is True:
#    os.remove(fastq_file)

#if sys.argv[1].endswith('.gz'):
#    p= subprocess.Popen('gzip %s' %(fastq_file), shell= True)
#    p.wait()

sys.exit()