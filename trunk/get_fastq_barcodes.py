#!/usr/local/bin/python

import operator
import sys

if len(sys.argv) == 1:
    sys.exit("""
Extracts the barcodes from a fastq file and prints them out in order of descending
frequency

USAGE
    get_fastq_barcodes.py <fastqfile>

OUTPUT

MEMO:
    Barcode is extracted from reads looking like this:
    @CRIRUN_795:3:1:2250:1032#ANNTNAA/1
""")

fh= open(sys.argv[1]) ## 'slx-4859.s_3_sequence.sanger.fq'
r= 0
bcode_dict= {}
for line in fh:
    if r % 4 == 0:
        bcode= line.strip()[-9:-2]
        bcode_dict[bcode]= bcode_dict.get(bcode, 0) + 1
        r += 1
    else:
        r += 1

sorted_x = sorted(bcode_dict.iteritems(), key=operator.itemgetter(1))
sorted_x.reverse()
n_reads= sum([x[1] for x in sorted_x])
print(sys.argv[1])
print('Total reads:    %s' %(n_reads))
print('Total barcodes: %s\n' %(len(sorted_x)))
for line in sorted_x:
    outline= '\t'.join([str(line[0]), str(line[1]), str(round(100 * (float(line[1])/n_reads), 2))])
    print(outline)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
