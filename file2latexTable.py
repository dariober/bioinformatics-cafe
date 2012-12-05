#!/usr/bin/env python

import sys
import subprocess
import tempfile
import os

docstring="""
DESCRIPTION

    Convert file tab separated, with header to latex table using R/xtable
    Latex code sent to stdout
USAGE
    file2latexTable.py <myfile.txt>
"""

if len(sys.argv) == 1 or sys.argv[1] in ['-h', '--help', '-help']:
    print(docstring)
    sys.exit()

tmpdir= tempfile.mkdtemp(suffix= 'file2latexTable')

texout= os.path.join(tmpdir, 'file2latexTable.tex')
cmd= '''
library(xtable)
tab<- read.table('%s', header= TRUE, sep= '\t', stringsAsFactors= FALSE)
xtab<- xtable(tab, caption= '', label= '')
print(file= '%s', x= xtab, include.rownames= FALSE, digits= 4)
''' %(sys.argv[1], texout)

rscript= os.path.join(tmpdir, 'file2latexTable.R')
rin= open(rscript, 'w')
rin.write(cmd)
rin.close()

p= subprocess.Popen('R CMD BATCH %s' %(rscript), shell= True)
p.wait()

print('')
for line in open(texout):
    print(line.strip())
print('')
sys.exit()