#!/usr/bin/env python

import sys
import subprocess
import os

docstring= """
DESCRIPTION
Convert a SAM file to BAM, sort and index.
NB: Input *.sam is replaced by *.bam

This is what is executed:
-------------------------
samtools view -S -u %(sam)s | samtools sort - %(bname)s &&
samtools index %(bname)s.bam &&
rm %(sam)s

USAGE
sam2bam.py <input.sam>

REQUIRES
samtools on PATH
"""

if len(sys.argv) != 2 or sys.argv[1] in ['-h', '--help', '-help']:
    sys.exit(docstring)

sam= sys.argv[1]
if not sam.endswith('.sam'):
    sys.exit('Input %s does not have extension .sam' %(sam))

## Create output file names:
bname= os.path.splitext(sam)[0]
unsortedbam= bname + ''

## NB: could use -u option but it throws and error (bug)

cmd= """samtools view -S -b %(sam)s | samtools sort - %(bname)s &&
samtools index %(bname)s.bam &&
rm %(sam)s
""" %{'sam': sam, 'bname': bname}

print(cmd)
p= subprocess.Popen(cmd, stdout= subprocess.PIPE, stderr= subprocess.PIPE, shell= True)
p.wait()
print(p.stderr.read())

sys.exit()