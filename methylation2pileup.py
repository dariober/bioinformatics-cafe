#!/home/berald01/.local/bin/python

import argparse
import os
import sys

parser = argparse.ArgumentParser(description= """

DESCRIPTION:

--------------------------------------------------------------------------------
Bismark methylation extractor version v0.7.4
CRIRUN_36:1:1:17428:1672#NGATGT/1       +       bsseq_synthetic 95      Z
CRIRUN_36:1:1:17428:1672#NGATGT/1       +       bsseq_synthetic 81      Z
CRIRUN_36:1:1:17428:1672#NGATGT/1       +       bsseq_synthetic 69      Z
CRIRUN_36:1:1:17428:1672#NGATGT/1       +       bsseq_synthetic 54      Z
CRIRUN_36:1:1:17428:1672#NGATGT/1       +       bsseq_synthetic 38      Z
CRIRUN_36:1:1:16447:1693#NGATGT/1       +       bsseq_synthetic 95      Z
--------------------------------------------------------------------------------
                 
OUTPUT:

--------------------------------------------------------------------------------
seqname    position    position+1    name    strand    nbases     calls    
chr1       0           1             chr1_0  +         10         Z:8,z:1,x:1

EXAMPLES:


TODO:

""", formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('input',
                    type= str,
                    help=""" Input file to pileup, typically produced by methylation_extractor
                    
                    """)
                    
parser.add_argument('--skip', '-s',
                    type= int,
                    default= 1,
                    help="""Skip this many lines of input. default is one to skip
the header line produced by methylation_extractor
                    
                    """)
                    
args = parser.parse_args()

" ----------------------------------------------------------------------------- "

nskip= 0
fh= open(args.input)
while nskip < args.skip:
    fh.readline()
    nskip += 1

piledict= {}
## This dict will be in the form: {'chr1':
##                                    {'position': {'X': 10, 'x': 20, ...}
##                                    }
##                                }

## my_dict[some_value] = my_dict.get(some_value, 0) + 1

for line in fh:
    line= line.strip().split('\t')
    chrom= line[2]
    position= int(line[3])
    call= line[4]
    if not chrom in piledict:
        piledict[chrom]= {}
    if not position in piledict[chrom]:
        piledict[chrom][position]= {}
    if not call in piledict[chrom][position]:
        piledict[chrom][position][call]= 1
    else:
        piledict[chrom][position][call] += 1

chroms= sorted(piledict.keys())

for c in chroms:
    positions= sorted(piledict[c].keys())
    for p in positions:
        calls= sorted(piledict[c][p].keys())
        print('\t'.join([c, str(p), str(piledict[c][p])]))
fh.close()

sys.exit()
