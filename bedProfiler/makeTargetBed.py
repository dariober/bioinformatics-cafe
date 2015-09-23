#!/usr/bin/env python

import pybedtools
import argparse
import sys
import random

VERSION='0.0.1'

parser = argparse.ArgumentParser(description= """DESCRIPTION
Generate target bed file for which features the profile is to be created.
Input features extended left and right by given amount. Extended intervals are
divided in bins of given size are renamed to have the central bin as 0.

Example:
echo "chr1 1000 1010" | tr ' ' '\t' \
| makeTargetBed.py -i - -g genome.txt -s 3 -c

chr1	1002	1003	-3
chr1	1003	1004	-2
chr1	1004	1005	-1
chr1	1005	1006	0
chr1	1006	1007	1
chr1	1007	1008	2
chr1	1008	1009	3
""", formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('--inbed', '-i', required= True, help='''Input bed file to use as target. Use - for reading from stdin.''')
parser.add_argument('--genome', '-g', required= True, help='''Genome file giving the size of the chromosomes. Format: chrom<tab>size''')
parser.add_argument('--center', '-c', action= 'store_true', help='''If set, input regions are reduced to a single-base interval in the middle of the region.
                    Maybe useful if the target bed has peaky features or motifs.''')
parser.add_argument('--slop', '-s', type= int, default= 1000, help='''Extend each target interval by this many bases left and right''')
parser.add_argument('--binSize', '-b', type= int, default= 1, help='''Divide the extened intervals in bins of this size''')

parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)

args= parser.parse_args()

def center(feature):
    """ Return interval feature reduced to a 1bp interval in the middle.
    """
    mid= round((feature.end - feature.start)/2)
    feature.start= feature.start + mid
    feature.end= feature.start + 1
    return feature

def renameBinToCenter(feature, slop, binSize):
    """Rename the bins assigend by windowMaker to subtract the slop size. In this
    way if the intervals are centered the mid point gets bin name 0 and left and
    right bins have -ve and +ve numbers, respectively.
    """
    feature.name= str(int(feature.name) - int(round(slop / binSize))-1)
    return feature

def shuffleTies(feature):
    """Shuffle the collapsed elements produced by groupBy -o collapse. Useful to
    pick one bin at random from overlapping ones.
    """
    ties= feature.name.split(',')
    random.seed(feature.name) # Random but same results from same input.
    feature.name= random.choice(ties)
    return feature

# ==============================================================================

if args.inbed == '-':
    inbed= sys.stdin
else:
    inbed= args.inbed

if args.center:
    bed= pybedtools.BedTool(inbed).each(center)
else:
    bed= pybedtools.BedTool(inbed)
    
slopBed= bed.slop(g= args.genome, b= args.slop)
windowBed= pybedtools.BedTool() \
    .window_maker(b= slopBed.fn, w= args.binSize, i= 'winnum') \
    .each(renameBinToCenter, args.slop, args.binSize) \
    .sort() \
    .groupby(g= [1,2,3], c= 4, o= ['collapse']) \
    .each(shuffleTies)

for line in windowBed:
    sys.stdout.write(str(line))

sys.exit()