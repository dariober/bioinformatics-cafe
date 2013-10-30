#!/usr/bin/env python

import argparse
import sys
import gzip

parser = argparse.ArgumentParser(description= """

DESCRIPTION
    Merge (group) lines in bed files in sliding windows with optional step.
    Merging is done by *number* of lines not distance. 
    
EXAMPLE
    cat in.bed 
1     chr1	3204562	3207049
2     chr1	3411782	3411982
3     chr1	3660632	3661579
4     chr1	4280926	4283093
5     chr1	4333587	4340172
6     chr1	4341990	4342162
7     chr1	4341990	4342162
8     chr1	4342282	4342918
9     chr1	4342282	4342918
10    chr1	4350280	4350395
    ...
    mergeBedLines.py -i in.bed -n 5
    chr1    3204562 4340172 5       chr1_3204562_4340172    .
    chr1    4341990 4350395 5       chr1_4341990_4350395    .
    ...
    mergeBedLines.py -i in.bed -n 5 -s 2
    chr1    3204562 4340172 5       chr1_3204562_4340172    .
    chr1    4280926 4342918 5       chr1_4280926_4342918    .
    chr1    4341990 4399322 5       chr1_4341990_4399322    .
    ...

NOTE:
    There is no check wheather the file is sorted by position!



""", formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('--bed', '-i',
                   required= True,
                   help='''Input bed file with lines to group''')

parser.add_argument('--nlines', '-n',
                   required= True,
                   type= int,
                   help='''Number of lines to include in each group.
                   ''')

parser.add_argument('--step', '-s',
                   type= int,
                   required= False,
                   default= 0,
                   help='''Step size by these many these many lines. Default 0: 
windows are not ovconsecutive without overlap.
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
        step= args.nlines
    else:
        step= args.step
    # --------------------------------------------------------------------------

    lineGroup= []  
    curChrom= None
    for line in bed:
        line= line.strip().split('\t')
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
    
    if lineGroup != []:
        print(mergeWindow(lineGroup))
        
if __name__ == '__main__':
    main()
    sys.exit()
