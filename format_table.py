#!/usr/bin/env python

import sys
import argparse

parser = argparse.ArgumentParser(description= """

    Reformat a space or tab separated file to have each column aligned by inserting the
    appropriate number of spaces
    
    E.g. from
    fc.enriched.bed 314660500       mc.enriched.bed 49501300        18811600        5.98    295848900       94.02
    hmc.enriched.bed        319880500       fc.enriched.bed 314660500       167726500       52.43   152154000       47.57
    
    TO:
    fc.enriched.bed   314660500  mc.enriched.bed 49501300   18811600   5.98   295848900  94.02
    hmc.enriched.bed  319880500  fc.enriched.bed 314660500  167726500  52.43  152154000  47.57
    
    USAGE
        format_table.py <table-file>
    <table-file>: File to format. Use - to read from stdin (e.g. less myfile | format_table.py -)
    Minimum number of spaces to separate the columns. Default: 4
    
    Read ten lines starting from line 3:    
        format_table.py -m 10 -N 3 <table-file>
    
        
    TODO:
        - Read CSV files: Now a column is split at every space (tab or space, one or more times).
        - Allow custom separator (pass to line.split())
        - Allow rugged files.


""", formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('input',
                    type= str,
                    help="""File to format or - to read from stdin.
                    """)

parser.add_argument('-n', '--nsep',
                    type= int,
                    default= 4,
                    help="""Separate columns by at least these many blank spaces.
    Default 4.
                    """)

parser.add_argument('-s', '--sep',
                    type= str,
                    default= '\t',
                    help="""Columns in input are separated by this string.
    Default is tab '\t'. NB: On bash shell a tab is inserted with Ctrl-V TAB.
    Typing '\t' will be interpreted literally.
                    """)

parser.add_argument('-m', '--maxlines',
                    type= int,
                    default= 0,
                    help="""Maximum number of lines to read.
                    An integer <= 0 will read the entire file (in memory)
    Default is 0 (read entire file).
                    """)

parser.add_argument('-N', '--nstart',
                    type= int,
                    default= 1,
                    help="""Start reading the file from this line number. 1-based
    Meaning that --nstart 1 will read from the beginning (default)
                    """)

args = parser.parse_args()

# -----------------------------------------------------------------------------
if args.input == '-':
    fh= sys.stdin
else:
    fh= open(args.input)

spacer= ' ' * args.nsep 
# -----------------------------------------------------------------------------

## Go to the requested line number
n= 1
while n < args.nstart:
    fh.readline()
    n += 1

## Prepare list of lists from the table-file: 
table= []
maxcols= 0
n= 0
for line in fh:
    n += 1
    line= line.strip().split(args.sep)
    if len(line) > maxcols:
        maxcols= len(line)
    table.append(line)
    if (args.maxlines > 0) and (n >= args.maxlines):
        break

## Get max column size for each column:
colwidths= [0] * maxcols
for i in range(0, maxcols):
    for line in table:
        cell= line[i]
        if len(cell) > colwidths[i]:
            colwidths[i]= len(cell)

## Print put the with the correct spacing
for line in table:
    for i in range(0, len(line)):
        colOffset= colwidths[i]
        cell= line[i]
        nspace= colOffset - len(cell)
        line[i]= cell + ' '*nspace
    print(spacer.join(line))
fh.close()
sys.exit()
