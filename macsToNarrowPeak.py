#!/usr/bin/env python

import sys
import math

docstring="""DESCRIPTION
    Convert the peak file *.xls from MACS1.4 to narrowPeak format (see https://sites.google.com/site/anshulkundaje/projects/idr)
    NB: - Peak name is generated as <chr_start_end> therefire it will differ from *.bed file
        - Pvalue -10*log10(pvalue) is converted to log10(pvalue)
USAGE
    macsToNarrowPeak.py <macs.peak.xls> >  <macs.peak.narrowPeak>
"""

if len(sys.argv) != 2 or sys.argv[1] in ['-h', '--help']:
    print(docstring)
    sys.exit()

def macsToLog10(x):
    """Convert -10*log10(pvalue) to -log10(pavlue).
    x= float. Typically 7th column of *.xls file
    Return: float
    """
    log10x= x/10
    return(log10x)

xls= open(sys.argv[1])
header= False
for line in xls:
    line= line.strip()
    if line.startswith('#') or line == '':
        continue
    if not header:
        header= True
        continue
    line= line.split('\t')
    npline= line[0:3]
    name= '_'.join(line[0:3])
    npline.append(name)
    npline.append('0')
    npline.append('+')
    npline.append(line[7]) ## fold_enrich
    npline.append( str(macsToLog10((float(line[6])))) ) ## -log10(pvalue)
    npline.append(str( -math.log10( (float(line[8])+10**-30)/100) )) ## -log10(FDR)
    npline.append('0') ## Point source
    npline= '\t'.join(npline)
    print(npline)
xls.close()
sys.exit()