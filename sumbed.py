#!/usr/local/bin/python

import sys
import subprocess
import os

if len(sys.argv) != 2:
    sys.exit("""

Return the sum of the regions spanned in a bed file after having merged the
overlapping features.

Use - to read file from stdin

USAGE

sumbed.py <myfile.bed>
cat myfile.bed | sumbed.py -

E.g.
chr1 1 10
chr1 20 25

sumbed.py >>> (10-1) + (25-20)= 14

Note: Requires mergeBed and input file must fit in memory

""")


fin= sys.argv[1]

if fin == '-':
    " Read bed ffrom stdout "
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
    
ftmp= sys.argv[1] + '.sumbed.bed'  ## Merged input bed
cmd= 'mergeBed -i %s > %s' %(fin, ftmp)
p= subprocess.Popen(cmd, shell= True)
p.wait()
merged= open(ftmp).readlines()

sumbed= 0
n= 0
for line in merged:
    line= line.rstrip('\n\r').split('\t')
    sumbed += int(line[2]) - int(line[1])
    n += 1
print('Size spanned (bp):                %s' %(sumbed))
print('Number of features after merging: %s' %(n))

os.remove(ftmp)
try:
    os.remove(fouttmp)
except:
    pass
sys.exit()