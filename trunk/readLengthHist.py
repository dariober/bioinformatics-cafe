#!/usr/bin/env python

import gzip
import sys

docstring= """DESCRIPTION
Histogram of read length in fastq file. 
Use - to read from stdin. Input fastq can be gzipped, but it's slower.

Output to stdout with tab separated coulumns, starting from shortest read length
found. Header is included:
<length> <count> <pct> <cum_count> <cum_pct> <hist>

Last column is a rudimentary graphical representation of the percentage reads
in each length.

USAGE
readLengthHist.py <fastq>
or
zcat reads.fq.gz | readLengthHist.py -
"""

if len(sys.argv) != 2 or sys.argv[1] in ['-h', '--help']:
    print(docstring)
    sys.exit()

if sys.argv[1] == '-':
    fin= sys.stdin
elif sys.argv[1].endswith('.gz'):
    fin= gzip.open(sys.argv[1])
else:
    fin= open(sys.argv[1])
    
read_len_hist= {}
n= -1
tot= 0
seqmass= 0
for line in fin:
    if n % 4 == 0:
        l= len(line.strip())
        seqmass += l
        read_len_hist[l]= read_len_hist.get(l, 0) + 1
        tot+=1
    n += 1

len_range= range(min(read_len_hist.keys()), max(read_len_hist.keys())+1)

header= '\t'.join(['length', 'count', 'pct', 'cum_count', 'cum_pct', 'hist'])
print(header)
cumcnt= 0
for i in len_range:
    if read_len_hist.has_key(i):
        cnt= read_len_hist[i]
    else:
        cnt= 0
    pct= round(100*(float(cnt)/tot), 2)
    cumcnt += cnt
    cumpct= round(100*(float(cumcnt)/tot), 2)
    bars= '|' * int(round(pct))
    bars_left= ' ' * (100 - int(round(pct)) - 1) + '+' 
    print(str(i) + '\t' + str(cnt) + '\t' + str(pct) + '\t' + str(cumcnt) + '\t' + str(cumpct) + '\t' + bars + bars_left)
print('## Sequence mass (bp)\t%s' %(seqmass))
