#!/home/berald01/.local/bin/python

import pysam
import sys
import os
import argparse
import re

parser = argparse.ArgumentParser(description= """

DESCRIPTION:
    Applies some fixes to a bam file.
    Current fix(es) are:
    - Replace missing sequence quality with '!' * length-of-read
      This is required by RNA-SeQC_v1.1.5.jar/GATK but note that it is allowed
      by the bam specifications to have missing quality.

EXAMPLES:

TODO:

""", formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('input',
                    type= str,
                    help="""Input bamfile to fix
                    """)

parser.add_argument('output',
                    type= str,
                    default= None,
                    nargs='?',
                    help="""Output bamfile. Default is <input-name>.fix.bam.
                    """)

args = parser.parse_args()

# -----------------------------------------------------------------------------

bamfile= pysam.Samfile(args.input, "rb" )
if args.output is None:
    outbam= pysam.Samfile(re.sub('\.bam$', '.fix.bam', args.input), "wb", template=bamfile)
else:
    outbam= pysam.Samfile(args.output, "wb", template=bamfile )

# -----------------------------------------------------------------------------

for AlignedRead in bamfile:
    if AlignedRead.qual is None:
        AlignedRead.qual= '!'*AlignedRead.rlen
    outbam.write(AlignedRead)

bamfile.close()
outbam.close()
sys.exit()