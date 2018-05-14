#!/usr/bin/env python3

import argparse
import gzip
import os
import sys

from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL) 

VERSION= '0.1.0'
thisprog= '%s %s' %(os.path.basename(__file__), VERSION)

parser = argparse.ArgumentParser(description= """
DESCRIPTION
Interleave one or more pairs of fastq files. Input files may be compressed in
which case they must have extension '.gz' or '.bz2'.

This program is not fast, on gzip input is probably 5-10x slower than a bash
script using Unix command (for example
https://github.com/dariober/bioinformatics-cafe/tree/master/catInterleaveFastq).
However, this program performs some sanity check and if the output is streamed
through an aligner the slow down should be neglegible.

USAGE
catInterleaveFastq.py -1 L0_R1.fq.gz L1_R1.fq.gz ... -2 L0_R2.fq.gz L1_R2.fq.gz ...
""", formatter_class= argparse.RawTextHelpFormatter, prog= os.path.basename(__file__))

parser.add_argument('--first', '-1',
                   help= '''First-in-pair fastq files                   
                   ''',
                   nargs= '+',
                   required= True)

parser.add_argument('--second', '-2',
                   help= '''Second-in-pair fastq files                   
                   ''',
                   nargs= '+',
                   required= True)

parser.add_argument('--version', '-v', action='version', version='%(prog)s ' + VERSION)

args= parser.parse_args()

if len(args.first + args.second) != len(set(args.first + args.second)):
    sys.stderr.write('\nDuplicate file names found in input\n')
    sys.exit(1)
    
if len(args.first) != len(args.second):
    sys.stderr.write('\nDifferent number of first-in-pair and second-in-pair fastq files\n')
    sys.exit(1)

def xopen(x):
    """Return file handle to file x. x may be compressed
    """
    if x.lower().endswith('.gz'):
        fh= gzip.open(x, 'rt')
    elif x.lower().endswith('.bz2'):
        fh= bz2.open(x, 'rt')
    else:
        fh= open(x)
    return fh

def next_fq_read(fh):
    """Read four line (i.e. a fastq record) from file handle fh.
    """
    rname= fh.readline()
    if rname == '':
        return ''
    assert rname.startswith('@')

    seq= fh.readline()
    cmt= fh.readline()
    assert cmt.startswith('+')
    
    qual= fh.readline()
   
    assert len(seq) == len(qual)

    rec= [rname, seq, cmt, qual]
    for x in rec:
        assert x != ''
    return rec

def clean_read_name(x):
    """Get the name of the read from 'x' where x is the read name from the
    fastq and can contain things like spaces or /1 /2 prefixes which are not
    part of the read name.
    """
    if x[0] == '@':
        x= x[1:]
    if x.endswith('/1') or x.endswith('/2'):
        x= x[0:len(x)-2]
    if ' ' in x:
        x= x[0:x.index(' ')]
    if '\t' in x:
        x= x[0:x.index('\t')]
    assert x != ''
    return x

for i in range(0, len(args.first)):
    first= xopen(args.first[i])
    second= xopen(args.second[i])
    while True:
        r1= next_fq_read(first)
        r2= next_fq_read(second)
        if r1 == '' and r2 == '':
            first.close()
            second.close()
            break
        if (r1 == '' and r2 != '') or (r1 != '' and r2 == ''):
            sys.stderr.write('\nUnequal number of reads in pair %s and %s' %(args.first[i], args.second[i]))
            sys.exit(1)
        if clean_read_name(r1[0]) != clean_read_name(r2[0]):
            sys.stderr.write('\nDifferent read names in interleaved input. Got records %s and %s\n' % (str(r1), str(r2)))
            sys.exit(1)
        sys.stdout.write(''.join(r1))
        sys.stdout.write(''.join(r2))
sys.exit()
