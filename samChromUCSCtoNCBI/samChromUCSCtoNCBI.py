#!/usr/bin/env python

import sys
import re
import os

VERSION='0.1.0'

docstring= """
DESCRIPTION
Convert chrom names in sam input from UCSC nomenclature (chr1, chrX, chrY, ..., chrM)
to NCBI nomenclature (1, X, ..., M). Fails if chrom names other than the full
chromosomes are found. E.g. "random", "unknown" contigs are not accepted.

USAGE
Takes no arguments: read sam with header from stdin and write sam to stdout
Example:
samtools view -h aln.bam | samChromUCSCtoNCBI.py | samtools view -S - > aln.ncbi.sam

Version %s
""" %(VERSION)

if len(sys.argv) != 1:
    print docstring
    sys.exit()

def ucsc2ncbi(chrom):
    """Convert chrom from ucsc to ncbi
    """
    if chrom in ['*', '=']:
        return chrom # For unmapped reads and PE reads mapped to same chrom
  
    if not chrom.startswith('chr'):
        sys.stderr.write('Chrom name does not start with "chr": %s\n' %(chrom))
        sys.exit()    

    newchrom= chrom[3:]
    if newchrom == 'M':
        newchrom= 'MT'
    elif newchrom in ['X', 'Y']:
        pass #
    else:
        int(newchrom) # Make sure name is an int!
    return(newchrom)

headerDone= False
addPGtag= True
for line in sys.stdin:
    line= line.strip().split('\t')
    if line[0].startswith('@'):
        if line[0].startswith('@SQ'):
            ucsc= re.sub('^SN:', '', line[1])
            ncbi= ucsc2ncbi(ucsc)
            line[1]= 'SN:' + ncbi
            headerDone= True
    else:
        if not headerDone:
            sys.stderr.write("Input sam doesn't have header")
            sys.exit()
        if addPGtag:
            addPGtag= False
            pgline= '\t'.join(['@PG', 'ID:'+os.path.basename(__file__), 'VN:'+VERSION])
            print pgline
        line[2]= ucsc2ncbi(line[2])
        line[6]= ucsc2ncbi(line[6])
    print '\t'.join(line)

