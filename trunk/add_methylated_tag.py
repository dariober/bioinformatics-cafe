#!/usr/bin/env python

import argparse
import sys
import os
import pysam

parser = argparse.ArgumentParser(description= """

DESCRIPTION:
    Output is SAM format to stdout.
    Add an additional tag YM with: <length of the read>-<number of cytosines in the read>-<number of cytosines methylated> 
    E.g. YM:Z:95-5-3
    
    This script extract the XM tag from each read to make it easy to add filters.

MEMO:

EXAMPLE

""", formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('--bam', '-b',
                   required= True,
                   help='''BAM file to open or - to read from stdin.
                   ''')

parser.add_argument('--ex_unmet', '-x',
                   action= 'store_true',
                   help='''Do not output reads completely methylated (completely unconverted).
If a == b in YM:Z:x-a-b do not print out.
                   ''')
def metString(alignedRead):
    """Parse the XM tag from bismark alignment to return a tuple as:
    (Methylation string, read length, no. of cytosines, no. cytosines methylated).
    alignedRead is a read object from pysam.
    Example output:
    ('H.hHh.h..hH.h...h...................................Z..........................................', 95, 10, 4)
    """
    metstring= [x[1] for x in alignedRead.tags if x[0] == 'XM'][0]
    cnt_M= sum(c.isupper() for c in metstring)
    cnt_m= sum(c.islower() for c in metstring)
    cnt_tot= cnt_M + cnt_m
    met_tuple= (metstring, len(metstring), cnt_tot, cnt_M)
    return(met_tuple)

def metString2Tag(met_tuple):
    """Convert metyhlation tuple from metString() to sam tag YM
    Output is tuple for pysam (YM, 'a-b-c')
    """
    ym= ('YM', '-'.join([str(x) for x in met_tuple[1:]]))
    return(ym)
    
if __name__ == '__main__':
    args = parser.parse_args()
    # main()
    ## Memo: if args.bam is '-' pysam will read from stdin.
    samfile = pysam.Samfile(args.bam, "rb" )
    outfile = pysam.Samfile("-", "wb", template = samfile)
    for alignedRead in samfile:
        methylation= metString(alignedRead) ## 
        ym= metString2Tag(methylation)
        if args.ex_unmet and (methylation[2] == methylation[3]):
            """Use methylation tuple to filter reads according to methylation status"""
            continue
        else:
            alignedRead.tags= alignedRead.tags + [ym]
            outfile.write(alignedRead)
    samfile.close()
    outfile.close()
    sys.exit()