#! /usr/bin/env python

import pybedtools as pyb
import sys

version= '0.1.0'
docstring= """Intersect two bed files and return the numbers for a venn diagram.

Usage
vennBedTwoWay.py a.bed b.bed
>>> c(only_a: 115841, only_b: 111798, both: 66130)

Version %s""" %(version)

def venn2(a, b):
    """Return the numbers for a two way venn diagram. a and b are the names of
    bed files.
    Return: dict of int {a, b, ab}
    """ 
    A= pyb.BedTool(a)
    B= pyb.BedTool(b)
    aonly= (A-B).count()
    bonly= (B-A).count()
    both= (B+A).count()
#    print('Only %s: %s' %(a, aonly))
#    print('Only %s: %s' %(b, bonly))
#    print('Both: %s' %(both))
    return({'only_a': aonly, 'only_b': bonly, 'both': both})

def vdict2R(vdict):
    """Print the dict produced by venn2 in an R-friendly format as
    named vector.
    """
    xv= """c(only_a: %(only_a)s, only_b: %(only_b)s, both: %(both)s)""" %(vdict)
    return(xv)
    
if __name__ == '__main__':
    if len(sys.argv) != 3 or sys.argv[1] in ['-h', '--help']:
        print(docstring)
        sys.exit(1)

    vdict= venn2(sys.argv[1], sys.argv[2])
    print( vdict2R(vdict) )
    
sys.exit()