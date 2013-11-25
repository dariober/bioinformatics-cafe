#!/usr/bin/env python

docstring= """
DESCRIPTION
    Returns the alpha and beta parameters for a beta distribution
    given the mean and standard deviation.
    Alpha and beta are called shape1 and shape2 in R.

    Memo: Beta distribution has mu and sd:
    mu= alpha / (alpha + beta)
    sd= sqrt( (a*b) / ((a+b)^2 * (a+b+1)) )
    
    For example usage of sympy:
    http://stackoverflow.com/questions/9440337/solving-systems-of-equations-with-sympy

USAGE & EXAMPLE:
    getBetaParams.py <Mean> <StdDev>

    ## Which parameters give a beta distribution with mean= 0.5 and sd= 0.2
    getBetaParams.py 0.5 0.2
    2.40
    0.60
"""
import sys
import sympy

if len(sys.argv) != 3 or sys.argv[1] in ('-h', '--help'):
    sys.exit(docstring)

mu= float(sys.argv[1])
std= float(sys.argv[2])

if mu <= 0 or std <= 0 or mu >= 1:
    sys.exit('mu and sd must be > 0. mu')

a= sympy.S('a')
b= sympy.S('b')

equations= [
    sympy.Eq(mu, a  / (a + b)),
    sympy.Eq(std, sympy.sqrt((a * b) / ((a+b)**2 * (a + b + 1)) ))
    ]
solutions= sympy.solve(equations)

for x in solutions:
    ## Discard a:0, b:0 solution
    if x[a] == 0 and x[b] == 0:
        continue
    elif x[a] <= 0 or x[b] <= 0:
        sys.exit('Invalid solutions: %s' %(solutions))
    else:
        solution= x
print(x[a])
print(x[b])
sys.exit()



