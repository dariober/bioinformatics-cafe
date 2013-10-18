#!/usr/bin/env python

import pysam
import sys

docstring= """DESCRIPTION
    Summarize methylation in reads aligned by bismark.

USAGE
    summarizeReadMethylation.py <in.bam> <min-calls>
    in.bam:
       Input bam file
    min-calls:
        Skip reads containing less than these many cytosines. Default 0
OUTPUT
    Histogrom to stdout.
    Stdout with tab separted columns: <% methylated> <no. reads> <% total reads in this bin>
"""

if len(sys.argv) < 2 or len(sys.argv) > 3 or sys.argv[1] in ('-h', '--help'):
    sys.exit(docstring)

if len(sys.argv) == 2:
    min_calls= 0
else:
    min_calls= int(sys.argv[2])

samfile= pysam.Samfile( sys.argv[1], "rb" )

pct_hist= {}
for i in range(0, 101):
    pct_hist[float(i)]= 0

n_tot= 0
for alnRead in samfile:
    n_tot += 1
#    if n_tot % 1000000 == 0:
#        sys.stderr.write('%s reads\n' %(n_tot))
    methylDict= {'Z': 0, 'z': 0, 'X': 0, 'x': 0, 'H': 0, 'h': 0, 'U': 0, 'u': 0, '.': 0}
    tags= alnRead.tags
    xm= [x[1] for x in tags if x[0] == 'XM']
    if xm == [] or len(xm) > 1:
        sys.exit("\nZero or more than one XM tags found. Read was:\n%s\n" %(alnRead))
    ## [('NM', 1), ('XX', '70'), ('XM', 'X.....XH....H......XH.........................................H.H.HH..'), ('XR', 'CT'), ('XG', 'GA')]
    xm= xm[0]
    for x in xm:
        methylDict[x] += 1
    M= methylDict['Z'] + methylDict['X'] + methylDict['H'] + methylDict['U']
    m= methylDict['z'] + methylDict['x'] + methylDict['h'] + methylDict['u']
    if (M+m) <= min_calls:
        continue
    else:
        pct_M= round( (float(M)/(M+m) ) * 100, 0)
        pct_hist[pct_M] += 1

tot_reads_scored= 0
for x in pct_hist:
    tot_reads_scored += pct_hist[x]

for i in sorted(pct_hist.keys()):
    pct_bin= round((pct_hist[i] / float(tot_reads_scored)) * 100, 2)
    print('\t'.join([str(i), str(pct_hist[i]), str(pct_bin)]))
sys.stderr.write('Reads analysed: %s\n' %(tot_reads_scored))
sys.stderr.write('Total reads: %s\n' %(n_tot))
samfile.close()