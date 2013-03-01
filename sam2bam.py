#!/usr/bin/env python

import sys
import subprocess
import os
import argparse

parser = argparse.ArgumentParser(description= """
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

""", formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('input',
                   help='''Input file as produced by samtools mpileup. Use - to read from stdin.

                   ''')

parser.add_argument('--noidx',
                   action= 'store_true',
                   help='''Do not index the sorted output bam.

                   ''')

parser.add_argument('--force', '-f',
                   action= 'store_true',
                   help='''Overwrite bam files if it exists. Default is to fail in such
case.

                   ''')

args= parser.parse_args()
# -----------------------------------------------------------------------------

sam= args.input
if not sam.endswith('.sam'):
    sys.exit('Input file "%s" does not have extension .sam' %(sam))

## Create output file names:
bname= os.path.splitext(sam)[0]
unsortedbam= bname + ''
bam= bname + '.bam'

if os.path.exists(bam) and not args.force:
    sys.exit('Bam file "%s" already exists. Use -f to overwrite' %(bam))

## NB: could use -u option but it throws and error (bug)

cmd_bam= """samtools view -S -b %(sam)s | samtools sort - %(bname)s &&""" %{'sam': sam, 'bname': bname}
cmd_idx= """samtools index %(bname)s.bam &&"""  %{'bname': bname}
cmd_rm= """rm %(sam)s""" %{'sam': sam}

if args.noidx:
    cmd= '\n'.join([cmd_bam, cmd_rm])
else:
    cmd= '\n'.join([cmd_bam, cmd_idx, cmd_rm])

print(cmd)
p= subprocess.Popen(cmd, stdout= subprocess.PIPE, stderr= subprocess.PIPE, shell= True)
p.wait()
print(p.stderr.read())

sys.exit()