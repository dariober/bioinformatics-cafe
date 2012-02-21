#!/usr/local/bin/python

import sys

print(sys.argv)

if len(sys.argv) < 3 or len(sys.argv) > 4:
    sys.exit("""
Divides each bed feature in n windows of equal size (after rounding).
E.g n=5:
chr 1 50 gene1
---->
chr 1  10 gene1 w1
chr 11 20 gene1 w2
chr 21 30 gene1 w3
chr 31 40 gene1 w4
chr 41 50 gene1 w5

USAGE
    bed_window.py <n windows> <inbed> <n skip> > stdout
    
    argv[1]: Number of windows (required)
    argv[2]: Input bed file (required)
    argv[3]: Number of lines to skip (default 0)
    
            """)

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
    lst= range(lst[0], lst[1]+1)
    division = len(lst) / float(n)
    full= [lst[int(round(division * i)): int(round(division * (i + 1)))] for i in xrange(n)]
    extremes= [[x[0], x[-1]] for x in full]
    return(extremes)

inbed= open(sys.argv[2])
n_windows= int(sys.argv[1])
try:
    nskip= int(sys.argv[3])
except:
    nskip= 0

n= 0
for line in inbed:
    if n < nskip:
        n += 1
        continue
    line= line.strip().split('\t')
    windows= partition([int(line[1]), int(line[2])], n_windows)
    for i in range(0, len(windows)):
        new_line= [line[0], windows[i][0], windows[i][1]] + line[1:] + [i+1]
        print('\t'.join([str(x) for x in new_line]))
    n += 1
sys.exit()