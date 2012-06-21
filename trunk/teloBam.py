#!/home/berald01/.local/bin/python

docstring= """
Read a bamfile file and reports the number of telomeric sequences found.
Search for: TTAGGG and CCCTAA

USAGE
    teloBam.py input.bam
OUTPUT
    Tab delimited row:
    - Input filename
    - Number of reads in input bam
    - Number of reads containing the telomeric motif
    - Total number of matches
 
"""

import sys
import pysam

if len(sys.argv) == 1 or sys.argv[1] == '-h' or sys.argv[1] == '--help':
    sys.exit(docstring)

fin= sys.argv[1] ## '/lustre/sblab/berald01/repository/bam/el016_hf2_9.bwa.hg18.bam'
SEQ= 'TTAGGG'
RCSEQ= 'CCCTAA' ## Reverse complement
n= 0
nreads= 0
nmatch= 0
samfile= pysam.Samfile(fin, 'rb')
for line in samfile:
    n += 1
    seq= line.seq
    fmatch= seq.count(SEQ)
    rmatch= seq.count(RCSEQ)
    if fmatch > 0 or rmatch > 0:
        nreads += 1
    nmatch += fmatch
    nmatch += rmatch
#    if n > 10000:
#        break
print('\t'.join([fin, str(n), str(nreads), str(nmatch)]))
samfile.close()
