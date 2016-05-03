#!/bin/env python

import sys
import argparse
import os

parser = argparse.ArgumentParser(description= """
DESCRIPTION
Read a sorted, TAB separated file and output non-duplicate lines.
By default duplicates are identified by the first record in each line.

Example:
sort -k1,1 infile.txt | getUniqueLines.py

See also
https://github.com/dariober/bioinformatics-cafe/tree/master/getUniqueLines
""", formatter_class= argparse.RawTextHelpFormatter, prog= os.path.basename(__file__))

parser.add_argument('--input', '-i',
                   required= False,
		   default= '-',
                   help='''Input file, appropriately sorted. Use - to read from stdin (default)
                   ''')

parser.add_argument('--col', '-c',
                   required= False,
		   type= int,
                   default= 1,
                   help='''Column index to use as unique key, default 1.
                   ''')
parser.add_argument('--version', action='version', version='%(prog)s 0.2.0')

args= parser.parse_args()

# ------------------------------------------------------------------------------

if args.col < 1:
    sys.stderr.write('\nColumn index must be > 0. Got %s\n\n' %(args.col))
    sys.exit(1)

SEP= '\t' # In/output separator
IDX= args.col - 1

if args.input == '-':
    fin= sys.stdin
else:
    fin= open(args.input)

stack= []
for line in fin:
    line= line.strip().split(SEP)
    if len(stack) == 0 or stack[-1][IDX] == line[IDX]:
        # Put reads with the same name in the stack
        stack.append(line)
    elif len(stack) == 1:
        # The current read name is different from the read name in the stack.
        # I.e. the read name in the stack is unique, output the read in the stack
        print SEP.join(stack[0])
        stack= [line]
    elif len(stack) >= 2:
        # Two reads with the same name: Discard both as you can't tell from which genome they come from
        stack= [line]
    else:
        sys.stderr.write('Unexpected case')
        sys.exit(1)

if len(stack) == 1:
    print SEP.join(stack[0])

fin.close()

sys.exit()
