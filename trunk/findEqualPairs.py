#!/usr/bin/env python

import sys
import argparse
import gzip
import itertools
import Levenshtein

parser = argparse.ArgumentParser(description= """
DESCRIPTION
    Find reads pairs with equal sequence in a pair of fastq files. The two fastq files
    are supposed to come from a paired-end run and therefore the order of the reads
    is the same in the files.
    Output is a tab separated columns with:
    <read name> <sequence mate 1> <sequence mate 2> <sequence mate 1> <edit distance>

USAGE:
    findEqualPairs.py -f <fastq-1> <fastq-2> 

""", formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('--fastq', '-f',
                   required= True,
                   nargs=2,
                   help='''Input fastq files, paired.
                   ''')

parser.add_argument('--distance', '-d',
                   required= False,
                   type= int,
                   default= 0,
                   help='''Levenshtein distance to consider the two
sequence equal. Default 0 (=identical).
                   ''')

parser.add_argument('--suffix', '-s',
                   required= False,
                   type= str,
                   choices= ['strip', '_', 'na'],
                   default= 'strip',
                   help='''Read names might have a suffix separated by space from the
rest of the name. Some aligners (bwa) strip this suffix, some others (bowtie/bismark)
merge it with the name using underscore. What do you want to do with it? Options are:
'strip': Remove it (like bwa, default)
'_': Merge it with '_' like bismark.
'na': Do nothing, just leave it a it is.
''')

args = parser.parse_args()

def cleanName(name, opt):
    """Edit the string name using the option opt"""
    pos= name.find(' ')
    if pos == -1:
        return(name)
    if opt == 'strip':
        name= name[0:pos]
    elif opt == '_':
        name= name.replace(' ', '_')
    else:
        sys.exit('Unexepcted option %s' %(opt))
    return(name)

fastq1= args.fastq[0]
fastq2= args.fastq[1]

if fastq1.endswith('.gz'):
    fq1= gzip.open(fastq1)
else:
    fq1= open(fastq1)

if fastq2.endswith('.gz'):
    fq2= gzip.open(fastq2)
else:
    fq2= open(fastq2)

i= 1
matchList= []
for r1, r2 in itertools.izip(fq1, fq2):
    if i == 1:
        name1= r1.rstrip().lstrip('@')
        name2= r2.rstrip().lstrip('@')
        i += 1
    elif i == 2:
        seq1= r1.strip()
        seq2= r2.strip()
        levDist= Levenshtein.distance(seq1, seq2)
        if levDist <= args.distance:
            if args.suffix != 'na':
                name1= cleanName(name1, args.suffix)
                ## name2= cleanName(name2, args.suffix)
            matchList.extend([name1, seq1, seq2, str(levDist)])
            print('\t'.join(matchList))
            matchList= []
        i += 1
    elif i == 3:
        i += 1
    elif i == 4:
        i= 1    
    else:
        sys.exit('Unexpected iteration')
    
fq1.close()
fq2.close()
sys.exit()
