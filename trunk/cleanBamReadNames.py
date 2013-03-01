#!/usr/bin/env python

import pysam
import sys
import subprocess
import tempfile
import argparse
import re

parser = argparse.ArgumentParser(description= """
DESCRIPTION
    Remove the /1 or /2 suffix from paired read names.
    
    bowtie2 and possibly other aligners append /1 and /2 to the read name first
    and second in pair. This causes problems to tools like bamUtils which expect
    paired read names to be identical.
    
EXAMPLE
    ## Input is BAM output BAM
    cleanBamReadNames.py -i test.bam -b > test.clean.bam
    
    ## Read SAM note -h option to samtools and -S
    samtools view -h test.bam | cleanBamReadNames.py -S -i -
    
    ## Use with bamUtils
    cleanBamReadNames.py -b -i test.bam | bam clipOverlap --in -.bam --out test.clip.bam --stats

""", formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('--input', '-i',
                   required= True,
                   help='''Input bam file to clean. Use - to read from stdin.

''')

parser.add_argument('--outbam', '-b',
                   required= False,
                   action= 'store_true',
                   help='''Write to stdout as BAM. Defualt is to write SAM.
''')

parser.add_argument('--isam', '-S',
                   required= False,
                   action= 'store_true',
                   help='''Input stream is SAM. Default is BAM.
                   ''')

args = parser.parse_args()


if args.input.endswith('.sam'):
    samfile = pysam.Samfile(args.input, "r")
elif args.input.endswith('.bam'):
    samfile = pysam.Samfile(args.input, "rb")
elif args.input == '-' and args.isam:
    samfile = pysam.Samfile("-", "r")
elif args.input == '-' and not args.isam:
    samfile = pysam.Samfile("-", "rb")
else:
    sys.exit('Extension of input file must be .sam or .bam')

if args.outbam:
    outfile = pysam.Samfile( "-", "wb", template = samfile)
else:
    outfile = pysam.Samfile( "-", "w", template = samfile)
    
for read in samfile:
    if read.qname[-2:] in ('/1', '/2'):
        read.qname= read.qname[0:-2]
        outfile.write(read)
    else:
        outfile.write(read)

samfile.close()
outfile.close()
sys.exit()
    
    
    
    