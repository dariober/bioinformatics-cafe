#!/usr/bin/env python

import pysam
import sys
import argparse
import matplotlib.pyplot as plt

parser= argparse.ArgumentParser()

argparse.ArgumentParser(description= """
DESCRIPTION:
    Produce a table of methylation bias for each position in a read aligned by bismark.

MEMO:

EXAMPLE

""", formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('--bam', '-b',
                   required= True,
                   help='''BAM file to open or - to read BAM from stdin.
                   ''')

parser.add_argument('--stop_after', '-s',
                   required= False,
                   default= -1,
                   type= int,
                   help='''Stop after reading this many reads. Default -1.
                   ''')

args= parser.parse_args()

def getXM(tags):
    """Get XM tag from tags list of tuples
    """
    xm= None
    for x in tags:
        if x[0] == 'XM':
            xm= x[1]
            break
    return(xm)

samfile = pysam.Samfile(args.bam)

def incrementHist(readHist, xm, qlen):
    """
    """
    for i in range(0, len(xm)):
        if xm[i] == 'Z':
            readHist[i][0] += 1
        elif xm[i] == 'z':
            readHist[i][1] += 1
        else:
            pass
        readHist[i][2] += 1
    return(readHist)

def addPct(readHist, mate):
    i= 1
    for x in readHist:
        if (x[0] + x[1]) == 0:
            pct_met= 'NA'
        else:
            pct_met= str(round(100 * (float(x[0]) / (x[0] + x[1])), 2))
        x.append(pct_met)
        x.append(i)
        x.append(mate)
        i += 1
        
def printer(readHist, header= True):
    """Print to stdout the methylation histogram
    """
    if header:
        print("\t".join(["cnt_met", "cnt_unmet", "cnt_total", "pct_met", "pos", "mate"]))
    for x in readHist:
        print("\t".join([str(z) for z in x]))
        
        
## One inner list for each position in the read. Inner list is [# Methylated, # Unmethylated, # Total]
readHist_R1= [[0, 0, 0]]
readHist_R2= [[0, 0, 0]] 
n= 0
for alignedread in samfile:
    if args.stop_after > 0 and args.stop_after <= n:
        break
    n += 1
    if n % 250000 == 0:
        sys.stderr.write("%s reads processed\n" %(n))

    xm= getXM(alignedread.tags)
    
    while(len(readHist_R1) < len(xm)):
        readHist_R1.append([0, 0, 0])

    while(len(readHist_R2) < len(xm)):
        readHist_R2.append([0, 0, 0])

    if alignedread.is_reverse:
        xm= xm[::-1]
    if not alignedread.is_paired or alignedread.is_read1: # not alignedread.is_reverse or :
        incrementHist(readHist_R1, xm)
    elif alignedread.is_read2: # not alignedread.is_reverse or :
        incrementHist(readHist_R2, xm)
    else:
        sys.exit("Unexpected mate or pairing information")

sys.stderr.write("%s reads processed\n" %(n))

r= 1
for readHist in [readHist_R1, readHist_R2]:
    if r == 1:
        header= True
    else:
        header= False
    addPct(readHist, mate= r)
    printer(readHist, header= header)
    r += 1

y= [x[3] for x in readHist_R1]
x= range(1, len(readHist_R1) + 1)
plt.plot(x, y)
plt.savefig("test.pdf")

samfile.close()
# Read SAM from stdin 