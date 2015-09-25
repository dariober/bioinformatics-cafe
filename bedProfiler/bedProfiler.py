#!/usr/bin/env python

#!/usr/bin/env python

import pybedtools
import argparse
import sys
import random
import subprocess

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
parser.add_argument('--nTargets', '-n', type= int, help='''Number of targets (e.g. CTCF sites) to use as normalizing factor. 
        If not given it's calculated by counting the number of '0' bins in input region.''')
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

def binListContains(xstr, xbin= '0'):
    """Return True if the string `xstr` contains the bin `xbin`. The list 
    of bins overlapping this interval must be in 4th field (name), comma
    separated as produced by `groupBy -o collapse`
    """
    bins= xstr.split(',')
    has_bin= xbin in bins
    return has_bin

if not args.nTargets:
    if args.regions.endswith('.gz'):
        cmd= 'gunzip -c '
    else:
        cmd= 'cat '
    cmd += '%s | cut -f 5 | grep -w "0" | wc -l' %(args.regions)
    p= subprocess.Popen(cmd, shell= True, stdout= subprocess.PIPE, stderr= subprocess.PIPE)
    p.wait()
    out, err = p.communicate()
    if p.returncode != 0:
        sys.stderr.write(err)
        sys.exit(1)
    nTargets= int(out.strip())
else:
    nTargets= args.nTargets

tab= pybedtools.BedTool(regions).intersect(b= scoreBedGraph, wa= True, wb= True, sorted= True) \
    .each(prepareForSortGroup) \
    .sort() \
    .groupby(g= [1], c= [1, 9], o= ['count', 'sum'])

# This is another hack: Resort to have bin (1st col) sorted as numeric
cmd= 'sort -k1,1n -s %s' %(tab.fn)
p= subprocess.Popen(cmd, shell= True, stdout= subprocess.PIPE, stderr= subprocess.PIPE)

print '\t'.join(HEADER)
for line in p.stdout:
    line= line.strip().split('\t')
    outline= featureToOutput(line, nTargets)
    outline= '\t'.join(str(x) for x in outline)
    print outline

sys.exit()
