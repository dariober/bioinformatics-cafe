#!/bin/env python

import sys

VERSION= '0.1.0'

docstring= """
DESCRIPTION
Read from stdin a (sorted) file and output non-duplicate lines.
Duplicates are identified by the first record in each line, tab separated.

EXAMPLE USAGE
sort -k1,1 f.sam rc.sam | getUniqueLines.py 

Version %s 
""" %(VERSION)

if len(sys.argv) != 1:
    print docstring
    sys.exit()

SEP= '\t' # In- out-put separator
IDX= 0 # Column index to discriminating filed. I.e. 0 to use first column

stack= []
for line in sys.stdin:
    line= line.strip().split(SEP)
    if len(stack) == 0 or stack[-1][IDX] == line[IDX]:
        # Put reads with the same name in the stack
        stack.append(line)
    elif len(stack) == 1:
        # The current read name is different from the read name in the stack.
        # I.e. the read name in the stack is unique, output the read in the stack
        print SEP.join(stack[0])
        stack= [line]
    elif len(stack) == 2:
        # Two reads with the same name: Discard both as you can't tell from which
        # genome they come from
        stack= [line]
    elif len(stack) > 2:
        sys.stderr.write('\nError: At most two reads should have the same name. Got\n%s\n\n' %(stack))
        sys.exit(1)
    else:
        sys.stderr.write('Unexpected case')
        sys.exit(1)

if len(stack) == 1:
    print SEP.join(stack[0])

sys.exit()
