#!/home/berald01/.local/bin/python

import sys
import argparse
import random

parser = argparse.ArgumentParser(description="""

DESCRIPTION
    Divides each bed feature in n windows of equal size (after rounding).
    E.g n=5:
    chr 1 50 gene1
    ---->
    chr 1  10 gene1 w1
    chr 11 20 gene1 w2
    chr 21 30 gene1 w3
    chr 31 40 gene1 w4
    chr 41 50 gene1 w5

    If a feature has a span lower than the number of windows, the windows will be
    of size 1bp and some (randomly chosen) windows will be duplicated to make up
    for the missing ones. 

OUTPUT:
    It's not true BED:
    1) chrom
    2) window start
    3) window end
    4+) Original input line starting from column 2 (i.e. excluding chrom)
    last) Window number
    
    """, formatter_class= argparse.RawTextHelpFormatter)


parser.add_argument('--nwinds', '-w',
                   type= int,
                   required= True,
                   help='''Number of windows.
                                   
                   ''')

parser.add_argument('bed',
                   type= str,
                   help='''Input bed file
                                   
                   ''')
parser.add_argument('--nskip', '-n',
                   type= int,
                   required= False,
                   default= 0,
                   help='''Number of lines to skip (e.g. 1 for skipping the header
line). Default 0.
                                   
                   ''')

parser.add_argument('--reverse', '-r',
                   action= 'store_true',
                   help='''By default windows are numbered 1:nwinds starting from
the leftmost (lowest) coordinate regardless of wheteher the bed feature is on the + or
- strand. With -r flag the features on the - strand get windows numbered nwinds:1
starting from leftmost coord. Use this option for profiling promoters or TSSs (so that
window 1 is the most distant from the TSS and window nwind is on the TSS).
                                   
                   ''')
args = parser.parse_args()
# -----------------------------------------------------------------------------
def sample_wr(population, k):
    "Chooses k random elements (with replacement) from a population"
    n = len(population)
    _random, _int = random.random, int  # speed hack 
    result = [None] * k
    for i in xrange(k):
        j = _int(_random() * n)
        result[i] = population[j]
    return(result)

def partition(lst, n): 
    """
    This fun divides a list in nearly equal chunks
    See here http://stackoverflow.com/questions/2659900/python-slicing-a-list-into-n-nearly-equal-length-partitions
    
    lst: List of len 2 with start and end of the seq (2nd and 3rd col of bed files). E.g. [0, 100]
    n: Number of partitions
    
    Output:
       List of lists where each inner list has the start and end of the chunk
    E.g.
    >>> partition([0,50], 5)
    [[0, 9], [10, 19], [20, 30], [31, 40], [41, 50]]    
    """
    lst[0]= lst[0]+1 ## This is because bed is 0-based but the algorithm below not (check though)
    lst= range(lst[0], lst[1]+1) ## Expand the whole range
    orin= n
    if len(lst) < n:
        n= len(lst)
    division = len(lst) / float(n)
    ## Put separators. Now you have list of lists
    full= [lst[int(round(division * i)): int(round(division * (i + 1)))] for i in xrange(n)]
    
    ## Get the start and end of each nested list.
    ## Take 1 from the start positions to return to bed being 0 based.
    ## Note that the end of one bin will become the start of the next.
    extremes= [[x[0]-1, x[-1]] for x in full]
    if len(lst) < orin:
        samples= sample_wr(extremes, orin-len(lst))
        extremes= sorted(extremes + samples)
    ## Increase bins by 1 if span is 0:
    for i in range(0, len(extremes)):
        b= extremes[i]
        if b[0] == b[1]:
            b[1]= b[1] + 1
            extremes[i]= b
    return(extremes)
# ------------------------------------------------------------------------------

inbed= open(args.bed)

n= 0
for line in inbed:
    if n < args.nskip:
        n += 1
        continue
    line= line.strip().split('\t')
    windows= partition([int(line[1]), int(line[2])], args.nwinds)
    for i in range(0, len(windows)):
        if args.reverse and line[5] == '-':
            j= (args.nwinds - 1) - i
        else:
            j= i
        new_line= [line[0], windows[i][0], windows[i][1]] + line[1:] + [j+1]
        print('\t'.join([str(x) for x in new_line]))
    n += 1
sys.exit()