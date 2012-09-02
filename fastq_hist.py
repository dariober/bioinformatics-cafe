#!/home/berald01/.local/bin/python

import argparse
import os
import sys
import gzip
import operator

parser = argparse.ArgumentParser(description= """
DESCRIPTION
    Reads sequences from fastq file and returns a histogram of their frequency
    
EXAMPLE
    ## Top 15 most represented sequences in a sample of 1000 reads, sampled every 10th from myfastq.fq:
    fastq_hist.py -t 15 -n 1000 -S 10 myfastq.fq

    ## Read from stdin:
    less myfastq.fq | fastq_hist.py -t 10 -

    ## Combine with fastx_trimmer to trim sequences:
    fastx_trimmer -Q33 -l 20 -m 20 -i myfastq.fq | fastq_hist.py - -t 10

SAMPLE OUTPUT

    ## file: -; N. reads read-in: 26466; N. in dictionary: 26466
    ATTAATGCGTAATTATTATT	21057	79.6%
    ATTAATGTGTAATTATTATT	2252	8.5%
    ATTAATGCGTAATTATTCTT	554	2.1%
    ATTAATGCGTTATTATTATT	90	0.3%
    ATTAATGCGTAATTATTGTT	81	0.3%
    ATTAATGCGTAGTTATTATT	70	0.3%
    ATTAGTGCGTAATTATTATT	66	0.2%


DEPENDS-ON:

DEPENDS-ON-ME:

TODO:

""", formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('infile',
                   type= str,
                   help='''Fastq file to summarize. Use - to read from STDIN.
Gzip'd files (ending in .gz) don't need to be uncompressed.
                   ''')
parser.add_argument('--nreads', '-n',
                   required= False,
                   default= -1,
                   type= int,
                   help='''Store up to this many reads in the dictionary.
Default to no limit.
                   ''')
parser.add_argument('--skip', '-s',
                   required= False,
                   default= 0,
                   type= int,
                   help='''Skip this many reads from the beginning of the file. 
                   ''')
parser.add_argument('--sample', '-S',
                   required= False,
                   default= 1,
                   type= int,
                   help='''Get a read every so often. (E.g. -S 10 will get a
read every 10th). Default is 1 (every read).
                   ''')
parser.add_argument('--top', '-t',
                   required= False,
                   default= -1,
                   type= int,
                   help='''For reporting: Print out only the top t reads.
Defualt no limit (-1)
                   ''')

args = parser.parse_args()

if args.infile == '-':
    fastq= sys.stdin
elif args.infile.endswith('.gz'):
    fastq= gzip.open(args.infile)
else:
    fastq= open(args.infile)

rdict= {}
n= 0 ## Number of reads currently read
m= 0 ## No. reads in dictionary
i= 2

for line in fastq:
    i += 1
    if i % 4 == 0:
        n += 1
        if n <= args.skip:
            continue
        if n % args.sample == 0:
            seq= line.strip()
            rdict[seq]= rdict.get(seq, 0) + 1
            m += 1
        if m >= args.nreads & args.nreads > 0:
            break

## Order dict by decreasing value
sorted_x = sorted(rdict.iteritems(), key=operator.itemgetter(1), reverse= True)
p= 0
print('\n## file: %s; N. reads: %s; N. in dictionary: %s; N. distinct reads: %s' %(args.infile, n, m, len(sorted_x)))
for x in sorted_x:
    p += 1
    vcount= x[1]
    perc= round(float(vcount) / m, 3) * 100
    print('\t'.join([x[0], str(vcount), str(perc) + '%']))
    if p >= args.top & args.top > 0:
        break

fastq.close()
sys.exit()