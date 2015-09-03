#!/usr/bin/env python

import sys
import subprocess
import argparse
import os
import tempfile
import shutil
import atexit

VERSION= '0.1.0'

parser = argparse.ArgumentParser(description= """
DESCRIPTION
Sort chromosomes in input bed file in the same order as in the header of
reference bam file. Chroms in bed but not in bam go last by position.
By default bed records are not further sorted within chromosomes. 

REQUIRES:
* samtools on PATH
* Unix tools: cut, sort

""", formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('--inbed', '-i', default= "-",
                   help='''Input bed file use - for stdin.
                   ''')

parser.add_argument('--bam', '-b', required= True,
                   help='''Reference bam file.
                   ''')

parser.add_argument('--sort', '-s', action= 'store_true',
                   help='''If set, sort by bed records by start and end within chromosomes.
Default is to leave them as in input.
                   ''')

parser.add_argument('--verbose', '-V', action= 'store_true',
                   help='''Print to stderr some messages about what's being done. Useful for debugging
                   ''')

parser.add_argument('--version', '-v', action='version', version='%(prog)s %(VERSION)s' %{'prog': '%(prog)s', 'VERSION': VERSION})

args = parser.parse_args()

# ------------------------------------------------------------------------------
# Some checks
# -----------
assert os.path.isfile(args.bam)
if args.inbed != "-":
    assert os.path.isfile(args.inbed)

# Get chrom names from bam
# ------------------------
cmd= "samtools view -H %s | grep -P '^@SQ' | sed 's/.*\tSN://' | sed 's/\t.*//'" %(args.bam)
if args.verbose:
    sys.stderr.write("Getting chromosome names and order:\n %s\n" %(cmd))
p= subprocess.Popen(cmd, shell= True, stdout= subprocess.PIPE, stderr= subprocess.PIPE)

chrom_names= p.stdout.readlines()
chrom_names= [x.strip() for x in chrom_names]
if len(chrom_names) == 0:
    sys.stderr.write("No chromsome found in input sam!")
if args.verbose:
    sys.stderr.write("Chromosomes found:\n %s\n" %(chrom_names))

# Dict of key/value as name/index
# -------------------------------
chrom_idx= {}
i= 0
for chrom in chrom_names:
    chrom_idx[chrom]= i
    i+=1
last_chrom_idx= i # For chroms in bed but in bam give an additional index which will place these chromosomes last.

if args.verbose:
    sys.stderr.write("Chromosomes (key) and indexes (value):\n %s\n" %(chrom_idx))

# Prepened to each bed record the chrom index
# -------------------------------------------
tmpdir= tempfile.mkdtemp(prefix= 'tmp_sortBedAsBam_')
atexit.register(shutil.rmtree, tmpdir)
tmpname= os.path.join(tmpdir, "inbed.txt")
tmpf= open(tmpname, 'w')
if args.verbose:
    sys.stderr.write("Writing tmp file to:\n %s\n" %(tmpname))

if args.inbed == "-":
    fin= sys.stdin
else:
    fin= open(args.inbed)
for line in fin:
    chrom= line.split('\t')[0]
    try:
        x_idx= chrom_idx[chrom]
    except KeyError:
        # This is for chroms in bed but not in bam
        chrom_idx[chrom]= last_chrom_idx
        x_idx= last_chrom_idx
        last_chrom_idx+=1
    outline= str(x_idx) + "\t" + line
    tmpf.write(outline)
fin.close()
tmpf.close()

# Sort file by chrom index and remove index column
if args.sort:
    cmd_sort= "sort -s -k1,1n -k3,3n -k4,4n "
else:
    cmd_sort= "sort -s -k1,1n "

cmd= cmd_sort + tmpname + " | cut -f2-"
if args.verbose:
    sys.stderr.write(cmd + "\n")

p= subprocess.Popen(cmd, shell= True, stdout= subprocess.PIPE, stderr= subprocess.PIPE)
for line in p.stdout:
    sys.stdout.write(line)

sys.exit(0)
