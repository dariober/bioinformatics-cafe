#!/usr/bin/env python

import sys
import os
import gzip

fq= sys.argv[1]
if fq.endswith('.gz'):
    fin= gzip.open(fq)
else:
    fin= open(fq)

trimer_dict= {}
i= -1
n= 0
for line in fin:
    if i % 4 == 0:
        n += 1
        trimer= line[0:3]
        if trimer in trimer_dict:
            trimer_dict[trimer] += 1
        else:
            trimer_dict[trimer]= 0
    i += 1
#    if i >= 1000:
#        break
trimer_list= []
for k in trimer_dict:
    trimer_list.append([trimer_dict[k], k])
trimer_list.sort(reverse= True)
trimer_list= ['\t'.join([str(x) for x in tup]) for tup in trimer_list]
trimer_list= [os.path.split(fq)[1]] + ['\t'.join([str(n), 'tot'])] + trimer_list
trimer_list= '\t'.join(trimer_list )

print(trimer_list)
fin.close()
