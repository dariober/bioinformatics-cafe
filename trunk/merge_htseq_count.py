#!/usr/bin/env python

import sys

docstring= """
DESCRIPTION
Sum counts in multiple file produced by htseq-count

USAGE
merge_htseq_count.py <file.1.htseq> <file.2.htseq> ...

Htseq file are tab separated with columns: <gene name> <count>.
They MUST have the same gene names in the same order. htseq-counts produces files
in this way.
"""

if (len(sys.argv) == 1) or (sys.argv[1] in ['-h', '--help']):
    print(docstring)
    sys.exit(1)

htseq= []
for f in sys.argv[1:]:
    sys.stderr.write('Opening "%s"\n' %(f))
    htseq.append(open(f))

n= 0
while True:
    first= True
    for f in htseq:
        line= f.readline().strip().split('\t')
        if line == ['']:
            sys.stderr.write("N lines: %s\n" %(n))
            for f in htseq:
                f.close()
            sys.exit()
        gene= line[0]
        cnt= int(line[1])
        if first:
            sumline= [gene, cnt]
            first= False
        elif sumline[0] != gene:
            sys.stderr.write("Unmatched gene names!\n")
            sys.exit(1)
        else:
            sumline[1] += cnt
    sumline[1]= str(sumline[1])
    print('\t'.join(sumline))
    n+=1

