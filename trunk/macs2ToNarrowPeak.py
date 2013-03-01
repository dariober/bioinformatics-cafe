#!/usr/bin/env python

import sys
import math

docstring="""DESCRIPTION
    Convert the peak file *.xls from MACSv2 to narrowPeak format (see https://sites.google.com/site/anshulkundaje/projects/idr)
    NB: - Peak name is generated as <chr_start_end> therefire it will differ from *.bed file
        - Pvalue -10*log10(pvalue) is converted to log10(pvalue)
USAGE
    macsToNarrowPeak.py <macs.peak.xls> >  <macs.peak.narrowPeak>

chr     start   end     length  abs_summit      pileup  -log10(pvalue)  fold_enrichment -log10(qvalue)  name
chr1    197     346     150     267     13.00000        15.77210        9.55000 -1.00000        el001_6.clean.hf2_peak_1
chr1    2803291 2803385 95      2803334 6.00000 5.07817 4.13396 -1.00000        el001_6.clean.hf2_peak_2
chr1    2957926 2958084 159     2958068 5.00000 4.06826 3.54340 -1.00000        el001_6.clean.hf2_peak_3
chr1    3316105 3316356 252     3316181 8.00000 7.26163 5.31509 -1.00000        el001_6.clean.hf2_peak_4
chr1    4104587 4104675 89      4104614 16.00000        8.15292 4.50550 -1.00000        el001_6.clean.hf2_peak_5
chr1    4831478 4831602 125     4831536 12.00000        8.50622 5.44712 -1.00000        el001_6.clean.hf2_peak_6
chr1    5120070 5120228 159     5120139 10.00000        6.58917 4.60910 -1.00000        el001_6.clean.hf2_peak_7
chr1    6013288 6013367 80      6013305 6.00000 5.07817 4.13396 -1.00000        el001_6.clean.hf2_peak_8
chr1    7086115 7086201 87      7086175 6.00000 5.07817 4.13396 -1.00000        el001_6.clean.hf2_peak_9
c
"""

if len(sys.argv) != 2 or sys.argv[1] in ['-h', '--help']:
    print(docstring)
    sys.exit()

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
    npline.append('.')
    npline.append(line[7]) ## fold_enrich
    npline.append(line[6]) ## -log10(pvalue)
    npline.append(line[8]) ## -log10(FDR)
    npline.append('-1') ## Point source
    npline= '\t'.join(npline)
    print(npline)
xls.close()
sys.exit()