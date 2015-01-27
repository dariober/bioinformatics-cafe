#! /usr/bin/env python

import pybedtools as pyb
import sys

version= '0.1.0'
docstring= """Intersect two bed files and return the numbers for a venn diagram.

Usage
vennBedTwoWay.py a.bed b.bed

Version %s""" %(version)

def venn2(a, b):
    """Return the numbers for a two way venn diagram. a and b are the names of
    bed files.
    Return: List of int [Only in a, Only in b, In both]
    """ 
    A= pyb.BedTool(a)
    B= pyb.BedTool(b)
    aonly= (A-B).count()
    bonly= (B-A).count()
    both= (B+A).count()
    print('Only %s: %s' %(a, aonly))
    print('Only %s: %s' %(b, bonly))
    print('Both: %s' %(both))
    return([aonly, bonly, both])

if len(sys.argv) < 3 or sys.argv[1] in ['-h', '--help']:
    print(docstring)
    sys.exit(1)

venn2(sys.argv[1], sys.argv[2])

sys.exit()