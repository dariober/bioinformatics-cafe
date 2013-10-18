#!/usr/bin/env python

import sys
import operator

bdict= {}
nlines= 0
for line in sys.stdin:
    nlines += 1
    line= line.split('\t')
    rname= line[0]
    bcode= rname[rname.find('#')+1:rname.find('#')+7]
    bdict[bcode]= bdict.get(bcode, 0) + 1
    if nlines % 1000000 == 0:
        sys.stderr.write(str(nlines) + '\n')

sorted_x = sorted(bdict.iteritems(), key=operator.itemgetter(1), reverse= True)

for line in sorted_x:
    print(line[0] + '\t' + str(line[1]))

#for bcode in sorted(bdict.keys(), reverse= True):
#    print(bcode + '\t' + str(bdict[bcode]))

sys.exit()
