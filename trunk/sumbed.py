#!/bin/env python

import sys
import subprocess
import os
import pybedtools

if len(sys.argv) != 2 or sys.argv[1] == '-h':
    sys.exit('''

Return the sum of the regions spanned in a bed file after having merged the
overlapping features.

Output is (tab separated):
    <input filename>
    <sum of regions *merged*>
    <number of features after merging>
    <sum of regions *w/o merging*>
    <number of features w/o merging (i.e wc)>

Use - to read file from stdin

USAGE

sumbed.py <myfile.bed>
cat myfile.bed | sumbed.py -
## Skip header line
tail -n+2 <myfile.bed> | sumbed.py -

E.g.
chr1 1 10
chr1 20 25

sumbed.py >>> (10-1) + (25-20)= 14

Note: Requires mergeBed and input file must fit in memory

''')
# ------------------------------------------------------------------------------
def sumbed(bedlist):
    """ Sum the features in a bed file to compute the overall span
    bedlist: list of list where each inner list is a line of a bed file (string).
    
    Returns a list of length 2: ['span bp', 'n lines']
    """
    span= 0
    n= 0
    for line in bedlist:
        span += int(line[2]) - int(line[1])
        n += 1
    return([span, n])
# ------------------------------------------------------------------------------
fin= sys.argv[1]

if fin == '-':
    " Read bed from stdout "
    fouttmp= 'sumbed.tmp.bed'
    fh= open(fouttmp, 'w')
    filein= sys.stdin
    while True:
        line= filein.readline()
        fh.write(line)
        if line == '':
            fh.close()
            break
    fin= fouttmp
    
## ftmp= sys.argv[1] + '.sumbed.bed'  ## Merged input bed
## cmd= 'mergeBed -i %s > %s' %(fin, ftmp)
## p= subprocess.Popen(cmd, shell= True)
## p.wait()
## bed= open(ftmp).readlines()

# ----------------------
## Stats for merged file
## ---------------------
bedin= pybedtools.BedTool(fin) ## Read-in input bed file
bedmerged= bedin.merge()## .merge() ## Merge it

bed= []
for line in bedmerged:
    bedline= []
    for x in line:
        bedline.append(x)
    bed.append(bedline)

sumsmerged= sumbed(bed)

## Stats for original file
bedOri= open(fin).readlines()
bedList= []
for line in bedOri:
    line= line.rstrip('\n\r').split('\t')
    bedList.append(line)

sums= sumbed(bedList)

print('\t'.join([fin, str(sumsmerged[0]), str(sumsmerged[1]), str(sums[0]), str(sums[1])]))

#print('Size spanned (bp):                %s' %(sumbed))
#print('Number of features after merging: %s' %(n))

#os.remove(ftmp)
#try:
#    os.remove(fouttmp)
#except:
#    pass
sys.exit()
