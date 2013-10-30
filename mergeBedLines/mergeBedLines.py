#!/usr/bin/env python

import argparse
import sys
import gzip

parser = argparse.ArgumentParser(description= """

DESCRIPTION
    Merge (group) lines in a bed file in sliding windows with optional step.
    Merging is done by *number* of lines not distance. 
    
EXAMPLE
    mergeBedLines.py -i sample.bed -n 10 -s 5

    chr1	0	38	chr1_0_38	10	.
    chr1	20	50	chr1_20_50	8	.
    chr2	100	122	chr2_100_122	6	.
    chr3	1000	1020	chr3_1000_1020	10	.
    chr3	1010	1022	chr3_1010_1022	6	.

NOTE:
    There is no check whether the file is sorted by position!
    You can check this with:
    sort -c -k1,1 -k2,2n -k3,3n sample.bed


""", formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('--bed', '-i',
                   required= True,
                   help='''Input bed file with lines to group, can be gzipped. Use - to read from stdin.''')

parser.add_argument('--nlines', '-n',
                   required= True,
                   type= int,
                   help='''Number of lines to include in each group.
                   ''')

parser.add_argument('--step', '-s',
                   type= int,
                   required= False,
                   default= 0,
                   help='''Step forward by these many lines after each window. Default 0: 
windows are consecutive without overlap.
                   ''')

# ------------------------------------------------------------------------------

def mergeWindow(lineList):
    """Return the merged list composed by the accumulated lines.
    lineList:
        List of list. Each inner list a line from the bed file
    Return:
        String
    """
    chrom= lineList[0][0]
    start= lineList[0][1]
    end= lineList[-1][2]
    name= chrom + '_' + start + '_' + end
    size= str(len(lineList))
    outline= '\t'.join([chrom, start, end, name, size, "."])
    return(outline)

def main():
    args = parser.parse_args()

    if args.bed == '-':
        bed= sys.stdin
    elif args.bed.endswith('.gz'):
        bed= gzip.open(args.bed)
    else:
        bed= open(args.bed)

    if args.nlines < args.step:
        sys.exit('Size of windows must be greater than step.')
    
    if args.step <= 0:
        "Step <= 0 means do not overlap. I.e. step == window size"
        step= args.nlines
    else:
        step= args.step
    # --------------------------------------------------------------------------

    lineGroup= []  
    curChrom= None
    n= 0
    for line in bed:
        line= line.strip().split('\t')[0:3]
        lineChrom= line[0]
        if curChrom is None:
            curChrom= lineChrom

        if len(lineGroup) >= args.nlines or curChrom != lineChrom:
            "Release lines once you have enough or the chrom ends"
            print(mergeWindow(lineGroup))
            lineGroup= lineGroup[step: ]
        
        lineGroup.append(line)

        if curChrom != lineChrom:
            curChrom= lineChrom
            lineGroup= [line]
        n += 1
        if n % 250000 == 0:
            sys.stderr.write(str(n) + ' processed\n')
        
    if lineGroup != []:
        print(mergeWindow(lineGroup))
    sys.stderr.write(str(n) + ' processed\n')
    
if __name__ == '__main__':
    main()
    sys.exit()
