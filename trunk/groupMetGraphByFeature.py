#!/usr/bin/env python

import sys
import subprocess
import gzip

docstring= """
DESCRIPTION
Summarizes the output in a bedgraph from methylation calling ("metGraph") by summing
the methylated and total counts by the features in a bed file

metGraph can be gzipped and it is expected to have columns:
<chrom> <start> <end> <pct met> <cnt met> <tot cnt> <strand>

Features bed is expected to have columns chrom, start, end (additional columns
ignored). Must be unzipped.

USAGE
groupBedgraphByFeature.py [-sorted] <features-bed> <metGraph>

-sorted:
    Use -sorted option in intersectBed. I.e. input files are sorted by position in the same way.
    There is no check for correct sorting.

OUTPUT FORMAT
Sorted by chrom, start, end:
<chrom> <feature start> <feature end> <pct met> <sum methylated counts> <sum tot counts> <strand or dot>
E.g.
chr1    135124  135563  0.85    17      20      .
chr1    713984  714547  0       0       112     .
chr1    762416  763445  0       0       203     .
"""

def execute(command):
    """ Execute shell command via subprocess and return output in real-time.
    See
    http://stackoverflow.com/questions/4417546/constantly-print-subprocess-output-while-process-is-running
    """
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    # Poll process for new output until finished
    while True:
        nextline = process.stdout.readline()
        if nextline == '' and process.poll() != None:
            break
        sys.stdout.write(nextline)
        sys.stdout.flush()

    output = process.communicate()[0]
    exitCode = process.returncode

    if (exitCode == 0):
        return output
    else:
        raise ProcessException(command, exitCode, output)

def checkMetGraph(metgraph):
    """Check bedgraph file has the expected format
    Expected:
    chr1    10496   10497   52.381  11      21      +
    chr1    10497   10498   0.0     0       2       -
    chr1    10524   10525   100.0   21      21      +
    chr1    10525   10526   100.0   21      21      -
    chr1    10541   10542   52.381  11      21      +
    """
    if metgraph.endswith('.gz'):
        bdg= gzip.open(metgraph)
    else:
        bdg= open(metgraph)
    n= 0
    for line in bdg:
        line= line.strip().split('\t')
        start= line[1]
        end= line[2]
        pct= line[3]
        cntM= line[4]
        tot=line[5]
        if len(line) != 7:
            sys.exit('Incorrect number of fields in %s' %(metgraph))
        if not all([start.isdigit(), end.isdigit(), cntM.isdigit(), tot.isdigit()]):
            sys.exit('Invalid line %s' %(metgraph))
        if int(cntM) < 0 or int(tot) < 0 or int(tot) < int(cntM):
            sys.exit("Invalid counts")
#       Check for pct deprecated as this filed is redundant and might be
#       replaced by some missing values in the future.
#        try:
#            float(pct)
#        except ValueError:
#            sys.exit("Invalid percentage format")
        if n > 100000:
            break
        n += 1
    bdg.close()

if '-sorted' in sys.argv:
    sorted= '-sorted'
    idx= sys.argv.index('-sorted') 
    del(sys.argv[idx])
else:
    sorted= ''

if len(sys.argv) != 3 or sys.argv[1] in ['-h', '--help']:
    print(docstring)
    sys.exit()
        
if sys.argv[1].endswith('.gz'):
    sys.exit("Feature file must be unzipped")

sys.stderr.write('Checking sample lines from bedgraph file...' + '\n')
checkMetGraph(sys.argv[2])

cmd= """
set -e
set -o pipefail

awk 'BEGIN{OFS= "\\t"} {print $1, $2, $3}' %s \\
| intersectBed -a - -b %s -wa -wb %s \\
| awk 'BEGIN {OFS= "\\t"} {print $1, $2, $3, $8, $9}' \\
| sort -k1,1 -k2,2n -k3,3n \\
| groupBy -i - -g 1,2,3 -c 4,5 -o sum,sum \
| awk 'BEGIN{OFS="\\t"} {print $1, $2, $3, $4/$5, $4, $5, "."}' \
""" %(sys.argv[1], sys.argv[2], sorted)

sys.stderr.write(cmd + '\n')
execute(cmd)

sys.exit()
