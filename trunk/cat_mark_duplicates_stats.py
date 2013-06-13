#!/usr/bin/env python

import glob
import sys
import os
import argparse

parser = argparse.ArgumentParser(description= """

DESCRIPTION

Concatenate duplicate stats from picard MarkDuplicates. Stats go to stdout, histogram to
file given via -H.

Both stats and hitogra, files will have filename left-bound.

USAGE
    cat_mark_duplicates_stats.py *.markDuplicates.txt -H hist.txt > concat.stats
    """, prog= 'cat_mark_duplicates_stats.py', formatter_class= argparse.RawDescriptionHelpFormatter) # , formatter_class= argparse.RawTextHelpFormatter

# -----------------------------------------------------------------------------

parser.add_argument('--version', action='version', version='%(prog)s 0.1.0')

parser.add_argument('--input', '-i',
                   required= True,
                   nargs= '+',
                   help='''Stat files to concatenate.
                   ''')

parser.add_argument('--hist', '-H',
                   required= False,
                    default= None,
                   help='''File with the concatenated *histograms*.
                   ''')


args = parser.parse_args()

def parseStats(x):
    """Extract line contaiining stats
    x:
        stat file name from MarkDuplicates
    return:
        Tuple with two lists: Header and stats
    """
    try:
        xstats= open(x)
    except IOError:
        sys.exit('I cannot open file %s' %x)
    xstats= xstats.readlines()
    header= None
    stats= None
    for i in range(0, len(xstats)):
        line= xstats[i]
        if line.startswith('## METRICS CLASS'):
            header= xstats[i+1].strip().split('\t')
            stats= xstats[i+2].strip().split('\t')
            break
    if not header or not stats:
        sys.exit('Cannot locate header and or stats line in file %s' %x)
    return((header, stats))

def parseHist(x):
    """Extract lines contaiining histogram
    x:
        stat file name from MarkDuplicates
    return:
        Tuple with two lists: Header and stats
    """
    try:
        xstats= open(x)
    except IOError:
        sys.exit('I cannot open file %s' %x)
    xstats= xstats.readlines()
    skip= None
    for i in range(0, len(xstats)):
        line= xstats[i]
        if line.startswith('## HISTOGRAM'):
            hist= xstats[i+2:]
    hist= [x.strip().split('\t') for x in hist if x.strip() != '']
    return(hist)

    
infiles= []
for f in args.input:
    infiles.extend(glob.glob(f))
infiles= sorted(list(set(infiles)))

if args.hist:
    hout= open(args.hist, 'w')
    hout.write('\t'.join(['FILENAME', 'BIN', 'VALUE']) + '\n')

header= False
for f in infiles:
    xf= parseStats(f)
    filename= xf[0]
    if not header:
        print('\t'.join(['FILENAME'] + filename))
        header= True
    if args.hist:
        xh= parseHist(f)
        for x in xh:
            x= [f] + x
            hout.write('\t'.join(x) + '\n')
        
    print('\t'.join([os.path.split(f)[1]] + xf[1]))

if args.hist:
    hout.close()
sys.exit()
