#!/usr/bin/env python

import pysam
import sys

VERSION= '0.1.0'

docstring= """Add tag with alignment length to input bam file.

USAGE
addAlignmentLengthToBam.py <in-bam> <tag name>

in bam
    Input bam or sam with header, use - to read from stdin.
tag name
    Name of tag to store length. Default 'XL'. Existing tags will be overwritten

Requires: pysam 0.8.3+

Version: %s
""" %(VERSION)

if len(sys.argv) < 2 or sys.argv[1] in ['-h', '--help']:
    print docstring
    sys.exit()

if sys.argv[1] in ['-v', '--version']:
    print VERSION
    sys.exit()

# ------------------------------------------------------------------------------

fname= sys.argv[1]

if len(sys.argv) == 2:
    tag= 'XL'
else:
    tag= sys.argv[2]
    
sam= pysam.AlignmentFile(fname)
out= pysam.AlignmentFile('-', 'wb', template= sam)
for line in sam:
    line.set_tag(tag, line.alen)
    out.write(line)

sam.close()
out.close()
sys.exit()