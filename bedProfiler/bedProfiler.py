#!/usr/bin/env python

#!/usr/bin/env python

import pybedtools
import argparse
import sys
import random

VERSION='0.0.1'

parser = argparse.ArgumentParser(description= """DESCRIPTION
Summarize features in bedGraph file around target features.

Example:
makeTargetBed.py -i - -g genome.txt -s 3 -c > bins.bed
bedProfiler.py -R bins.bed -s test_data/scores.bedGraph

TODO: Accept --regions from stdin
""", formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('--regions', '-R', required= True, help='''Input bed file with regions to plot. Typically produced by makeTargetBed.py''')
parser.add_argument('--scoreBedGraph', '-s', required= True, help='''Bedgraph file with scores to plot. - for stdin''')

parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)

args= parser.parse_args()

# ==============================================================================

# Output header
HEADER= ['bin', 'nHits', 'scoreSum', 'nTargets', 'scoreAvg']

regions= pybedtools.BedTool(args.regions)

if args.scoreBedGraph == '-':
    scoreBedGraph= pybedtools.BedTool(sys.stdin)
else:
    scoreBedGraph= pybedtools.BedTool(args.scoreBedGraph)

def prepareForSortGroup(feature):
    """This is a hack: We want to groupBy bin, which is in 4th column ("name").
    For this we need to sort by col 4, which is not possible with sortBed. So
    replace chrom with bin and sort by position.
    Really, you should avoid sortBed and use gnu sort instead.
    """
    feature.chrom= feature.name
    return feature

def featureToOutput(feature, nTargets):
    """Compute average of score by dividing the sum by the number of targets.
    Returns list with nTargets and average.
    feature:
        List with elements: [bin:str, nHits:int, avgSum:float]
    nTargets:
        Denominator for average, typically the number of features
    """
    avg= float(feature[2]) / nTargets
    outlst= [feature[0], feature[1], feature[2], nTargets, avg]
    return outlst
    
nTargets= regions.filter(lambda x: x.name == '0').count()

tab= pybedtools.BedTool(regions).intersect(b= scoreBedGraph, wa= True, wb= True, sorted= True) \
    .each(prepareForSortGroup) \
    .sort() \
    .groupby(g= [1], c= [1, 8], o= ['count', 'sum'])

print '\t'.join(HEADER)
for line in open(tab.fn):
    line= line.strip().split('\t')
    outline= featureToOutput(line, nTargets)
    outline= '\t'.join(str(x) for x in outline)
    print outline

#for feature in tab:
#    outline= featureToOutput(feature, nTargets)
#    print outline

"""
intersectBed -sorted -a targets.bed.gz -b $bdg -wa -wb \
| cut -f 4,8,9 \
| sort -k1,1n -s \
| awk -v OFS='\t' '{print \$1, 0, 0, \$2, \$3}' \
| groupBy -i - -g 1 -c 1,4,5 -o count,sum,sum \
| cut -f 1,2,3,4 \
| awk -v OFS='\t' -v nTargets=$nTargets '{print \$1, \$2, nTargets, \$3/nTargets, \$4/nTargets}' > $tab" > ${tab}.sh
"""