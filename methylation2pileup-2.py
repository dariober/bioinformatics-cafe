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
                    default= 0,
                    help="""Skip this many lines of input. default is zero.
                    
                    """)
                    
args = parser.parse_args()

" ----------------------------------------------------------------------------- "
def initline(chrom, pos, outputcalls):
    """ Initialize a line for output which will look like:
    line= ['chr1', '123', {'z': 0, 'Z': 0, 'x': 0, 'X': 0, 'h': 0, 'H': 0}]"""
    d= {}
    for x in outputcalls:
        d[x]= 0
    line= [chrom, pos, d]
    return(line)

def writeline(line, outputcalls):
    """ Print-out a line in the format produced by initline according to
    the order in outputcalls"""
    outline= [line[0], line[1]]
    for k in outputcalls:
        outline.append(str(line[2][k]))
    print('\t'.join(outline))
    

" ----------------------------------------------------------------------------- "
OUTPUTCALLS= ['z', 'Z', 'x', 'X', 'h', 'H'] ## Calls and their order to be sent to output

## Print header
print('\t'.join(['chr', 'pos'] + OUTPUTCALLS))

nskip= 0
fh= open(args.input)
while nskip < args.skip:
    fh.readline()
    nskip += 1

lastline= False
curline= fh.readline().strip().split('\t')
curchrom= curline[2]
curpos= curline[3]
call= curline[4]
outline= initline(curchrom, curpos, OUTPUTCALLS)
outline[2][call] += 1
while True:
    line= fh.readline().strip().split('\t')
    if line == ['']:
        break
    chrom= line[2]  
    pos= line[3]
    call= line[4]
    while chrom == curchrom and pos == curpos:
        outline[2][call] += 1
        line= fh.readline().strip().split('\t')
        if line == ['']:
            lastline= True
            break
        chrom= line[2]  
        pos= line[3]
        call= line[4]
    writeline(outline, OUTPUTCALLS)
    curchrom= chrom
    curpos= pos
    outline= initline(curchrom, curpos, OUTPUTCALLS)
    outline[2][call] += 1
if not lastline:
    writeline(outline, OUTPUTCALLS)

sys.exit()
