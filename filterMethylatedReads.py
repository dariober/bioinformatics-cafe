#!/usr/bin/env python

import sys
import pysam
import argparse


parser = argparse.ArgumentParser(description= """
DESCRIPTION
    Filter out reads from a bam file if they contain more than a threshold level of
    methylated reads in non-CpG context.
    
    Input is SAM or BAM file, typically produced by bismark, where the XM tag can
    be extracted (e.g: XM:Z:......z....xz...h..h.)

USAGE
    filterMethylatedReads.py -i <myreads.bam> -F 2 > filtered.bam

MEMO: Bismark codes for methylation are:
    ---------------------------------
    z   unmethylated C in CpG context
    Z   methylated C in CpG context
    x   unmethylated C in CHG context
    X   methylated C in CHG context
    h   unmethylated C in CHH context
    H   methylated C in CHH context
    ---------------------------------
    
""", formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('--input', '-i',
                   required= True,
                   help='''Input bam or sam file or - to read from stdin
                   ''')

parser.add_argument('--isam', '-S',
                   required= False,
                   action= 'store_true',
                   help='''Input stream is SAM. Default is BAM.
                   ''')

parser.add_argument('--outbam', '-b',
                   required= False,
                   action= 'store_true',
                   help='''Write to stdout as BAM. Defualt is to write SAM.
                   ''')

parser.add_argument('--filter', '-F',
                    required= True,
                    type= str,
                    help='''Filter reads containing methylated calls in non-CpG
context. If int, reads with more non-CpG calls than F will be filtered out.
If float between 0 and 1, reads with percentage of non-CpG calls greater than
F will be removed.
                   ''')

args = parser.parse_args()

# -----------------------------------------------------------------------------
def getTag(tagList, tag):
    """Return the index associated to to the tag from the list
    of tags in pysam.AlignedRead.tags
    E.g. getTag([('AS', -8), ('XN', 0), ('XM', 0)], 'AS') >>> 0
    """
    i= 0
    for t, v in tagList:
        if t == tag:
            return(i)
        i += 1
    return(None)

def countMethylation(xm):
    """Produce a dictionary of counts for each character in the XM tag
    countMethylation('......z....xz...h..h.') >>> {'z': 2, 'h': 2, 'X': 0, 'H': 0, 'x': 1, 'Z': 0}
    """
    metDict= {'z': xm.count('z'),
              'Z': xm.count('Z'),
              'x': xm.count('x'),
              'X': xm.count('X'),
              'h': xm.count('h'),
              'H': xm.count('H')}
    return(metDict)

def filterType(filter):
    """Return True if string 'filter' is int, False if filetr is a float
    between 0 and 1, None otherwise
filterType('1')   ## True
filterType('0')   ## True
filterType('1.0') ## False
filterType('0.0') ## False    
filterType('1.1') ## None
filterType('foo') ## None
filterType('-1') ## None
    """
    if filter.isdigit():
        if int(filter) < 0:
            return(None)
        else:
            return(True)
    else:
        try:
            filter= float(filter)
        except ValueError:
            return(None)
        if filter >= 0 and filter <= 1:
            return(False)
        else:
            return(None)
# -----------------------------------------------------------------------------

if args.input.endswith('.sam'):
    samfile = pysam.Samfile(args.input, "r")
elif args.input.endswith('.bam'):
    samfile = pysam.Samfile(args.input, "rb")
elif args.input == '-' and args.isam:
    samfile = pysam.Samfile("-", "r")
elif args.input == '-' and not args.isam:
    samfile = pysam.Samfile("-", "rb")
else:
    sys.exit('Extension of input file must be .sam or .bam')



if args.outbam:
    outfile = pysam.Samfile( "-", "wb", template = samfile)
else:
    outfile = pysam.Samfile( "-", "w", template = samfile)

isint= filterType(args.filter)
if isint is None:
    sys.exit('Invalid value for filter: "%s"' %(args.filter))

for read in samfile:
    xmidx= getTag(read.tags, 'XM')
    xm= read.tags[xmidx][1]
    cntDict= countMethylation(xm)
    cntXH= cntDict['X'] + cntDict['H']
    print(cntDict, cntXH)
    if isint:
        if cntXH <= int(args.filter):
            outfile.write(read)
    else:
        if (float(cntXH) / len(xm)) <= float(args.filter):
            outfile.write(read)
outfile.close()
samfile.close()

sys.exit()
    
    
    