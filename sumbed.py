#!/usr/bin/env python

import sys
import subprocess
import os
import pybedtools
import tempfile
import atexit

if len(sys.argv) != 2 or sys.argv[1] == '-h':
    sys.exit('''

Return the sum of the regions spanned in a bed file after having merged the
overlapping features.

Output is tab separated:
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

''')
# ------------------------------------------------------------------------------

class Sumbed:
    def __init__(self):
        self.filename= '-'
        self.spanMerged= 0
        self.nMerged= 0
        self.span= 0
        self.n= 0
    def incrementMerged(self, bedLine):
        """Increment counter for merged features.
        bedLine: A pybedtool.Interval
        """
        self.nMerged+= 1
        self.spanMerged+= bedLine.length
    def increment(self, bedLine):
        """Increment counter for original features.
        bedLine: A pybedtool.Interval
        """
        self.n+= 1
        self.span+= bedLine.length
    def toString(self):
        """Print object as string
        """
        tostr= self.filename + '\t' + \
               str(self.spanMerged) + '\t' + \
               str(self.nMerged) + '\t' + \
               str(self.span) + '\t' + \
               str(self.n)
        return(tostr)


# ------------------------------------------------------------------------------
infile= sys.argv[1]

tmp= tempfile.NamedTemporaryFile(prefix= 'tmp_sumbed_', suffix= '.bed', delete= False)
atexit.register(os.remove, tmp.name)

if infile == '-':
    " Read bed from stdout "
    fouttmp= tmp.name
    fh= open(fouttmp, 'w')
    filein= sys.stdin
    while True:
        line= filein.readline()
        fh.write(line)
        if line == '':
            fh.close()
            break
    fin= fouttmp
else:
    fin= infile

sumbed= Sumbed()
sumbed.filename= infile
    
# Stats for merged file
# ---------------------
bedin= pybedtools.BedTool(fin) 
bedmerged= bedin.sort().merge()

for line in bedmerged:
    sumbed.incrementMerged(line)

# For original file
# ----------------
bedOri= pybedtools.BedTool(fin)
for line in bedOri:
    sumbed.increment(line)

print sumbed.toString()

sys.exit()

