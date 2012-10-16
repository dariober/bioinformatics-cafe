#!/home/berald01/.local/bin/python

import argparse
import os
import sys
import decimal

parser = argparse.ArgumentParser(description= """

DESCRIPTION:

Convert output from bismark/methylation_extractor to produce a tabular/pileup file of
the methylation state at each C.
MEMO: Input is the output of methylation_extractor with reads sorted by chromosome
(sorting by position within chromosomes is not necessary a 'sort -k 3,3' should suffice.
Alternatively execute methylation_extractor on sorted sam file).

--------------------------------------------------------------------------------
CRIRUN_56:1:1:15122:1594#GCCAAT/1       -       chr1    103485851       x
CRIRUN_56:1:1:15122:1594#GCCAAT/1       -       chr8    103485851       X
CRIRUN_56:1:1:15122:1594#GCCAAT/1       -       chr8    103485855       X
...
CRIRUN_56:1:1:15122:1594#GCCAAT/1       -       chr8    103485868       H
CRIRUN_56:1:1:15122:1594#GCCAAT/1       -       chr8    103485868       X
CRIRUN_56:1:1:15122:1594#GCCAAT/1       -       chrX    103485889       h
CRIRUN_56:1:1:15122:1594#GCCAAT/1       -       chrM    103485868       H
CRIRUN_56:1:1:15122:1594#GCCAAT/1       -       chrM    103485878       h
...
CRIRUN_56:1:1:15122:1594#GCCAAT/1       -       chrM    103485889       h
--------------------------------------------------------------------------------
                 
OUTPUT (header line is included):

--------------------------------------------------------------------------------
chr	pos             z	Z	x	X	h	H
chr1	103485851	0	0	1	0	0	0
chr8	103485851	0	0	0	1	0	0
chr8	103485855	0	0	0	1	0	0
chr8	103485862	0	0	0	0	1	0
chr8	103485863	0	0	0	0	1	0
chr8	103485864	0	0	0	0	1	0
chr8	103485868	0	0	0	1	0	1
chrX	103485889	0	0	0	0	1	0
chrM	103485868	0	0	0	0	0	1
chrM	103485878	0	0	0	0	1	0
chrM	103485883	0	0	0	0	3	0
chrM	103485889	0	0	0	0	1	0
--------------------------------------------------------------------------------
EXAMPLES:


TODO:
- Check input is sorted by chromosome

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

def sumMethylation(calldict):
    """ Sum methylation calls in calldict and return percent metheylated
    calldict: {'z': 2, 'Z': 10, 'h': 0, 'H': 0, 'x': 0, 'X': 0}"""
    
    methylationSumDict= {'tot': 0, 'M':0 , 'pM': 0}
    
    for k in calldict:
        if k in ('Z', 'X', 'H'):
            methylationSumDict['M'] += calldict[k]
        methylationSumDict['tot'] += calldict[k]
    methylationSumDict['pM']= round(100*(float(methylationSumDict['M'])/(methylationSumDict['tot'])), 3) ## % methylated
    return(methylationSumDict)

def writeline(chrom, chrdict, outputcalls):
    """chrom= str for the name of the current chrom
    chrdict= dict of dicts with calls at each position: chrdict= {1:  {'z': 2, 'Z': 10, 'h': 0, 'H': 0, 'x': 0, 'X': 0}, 
                                                                  4:  {'z': 0, 'Z': 0, 'h': 0, 'H': 0, 'x': 3, 'X': 20}, 
                                                                  10: {...}
                                                                 }
    outputcalls: List of tags to represent each call from bismark: ['z', 'Z', 'x', 'X', 'h', 'H', ...] 
    """
    sortpos= sorted(chrdict.keys())
    for x in sortpos:
        callLine= chrdict[x]
        lineSum= sumMethylation(callLine) ## Add summary of meyhtaltion
        for k in lineSum:
            callLine[k]= lineSum[k]
        outline= [chrom, str(x)]
        for k in outputcalls:
            outline.append(str(callLine[k]))
        print('\t'.join(outline))
    
" ----------------------------------------------------------------------------- "

decimal.getcontext().prec = 2

OUTPUTCALLS= ['z', 'Z', 'x', 'X', 'h', 'H'] ## Calls and their order to be sent to output
outputcallsSum= OUTPUTCALLS + ['tot', 'M', 'pM']
## Print header
print('\t'.join(['chr', 'pos'] + outputcallsSum))

nskip= 0
fh= open(args.input)
while nskip < args.skip:
    fh.readline()
    nskip += 1

curline= fh.readline().strip().split('\t')
curchrom= curline[2]
curpos= int(curline[3])
chrdict= {curpos: {}}
chrdict[curpos]= {}
for x in OUTPUTCALLS:
    chrdict[curpos][x]= 0
call= curline[4]
chrdict[curpos][call] += 1
while True:
    line= fh.readline().strip().split('\t')
    if line == ['']:
        break
    chrom= line[2]  
    pos= int(line[3])
    call= line[4]
    while chrom == curchrom:
        if not pos in chrdict:
            chrdict[pos]= {}
            for x in OUTPUTCALLS:
                chrdict[pos][x]= 0
        chrdict[pos][call] += 1
        line= fh.readline().strip().split('\t')
        if line == ['']:
            break
        chrom= line[2]  
        pos= int(line[3])
        call= line[4]
    if line == ['']:
        writeline(curchrom, chrdict, outputcallsSum)
        break
    writeline(curchrom, chrdict, outputcallsSum)
    curchrom= chrom
    curpos= pos
    chrdict= {}
    chrdict[curpos]= {}
    for x in OUTPUTCALLS:
        chrdict[curpos][x]= 0
    call= line[4]
    chrdict[curpos][call] += 1
sys.exit()
