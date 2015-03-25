#! /usr/bin/env python

import pybedtools as pyb
import sys

version= '0.1.0'
docstring= """Intersect two bed files and return the numbers for a venn diagram.

Usage
vennBedTwoWay.py a.bed b.bed
>>> c(only_a: 115841, only_b: 111798, both: 66130)

Within R:
xv<- system('vennBedTwoWay.py a.bed b.bed', intern= TRUE)
ab<- eval(parse(text= xv))

See also https://pythonhosted.org/pybedtools/3-brief-examples.html

*Important note* The intersection count is not symmetrical! I.e. A x B != B x A
so the output here should be taken with care.

A.bed  ----    -------------        
B.bed           ---     ---   ----

There is only-A: 1; only-B: 1; both: 1 if A x B but 2 if B x A.

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
    both= (A+B).count()
    vdict= {'only_a': aonly, 'only_b': bonly, 'both': both}
    return(vdict)

def vdict2R(vdict):
    """Print the dict produced by venn2 in an R-friendly format as
    named vector.
    """
    xv= """c(only_a= %(only_a)s, only_b= %(only_b)s, both= %(both)s)""" %(vdict)
    return(xv)
    
if __name__ == '__main__':
    if len(sys.argv) != 3 or sys.argv[1] in ['-h', '--help']:
        print(docstring)
        sys.exit(1)

    vdict= venn2(sys.argv[1], sys.argv[2])
    print( vdict2R(vdict) )
    
sys.exit()
