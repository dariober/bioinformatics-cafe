#!/usr/bin/env python

import pysam
import sys
import os
import argparse

VERSION= '0.1.0'
PN= os.path.basename(__file__)

parser= argparse.ArgumentParser(formatter_class= argparse.RawTextHelpFormatter, description ="""
Add tags with alignment length, start and end to input bam file.
See http://pysam.readthedocs.org/en/latest/api.html#pysam.AlignedSegment for
details of tags. Requires: pysam 0.8.3+

For all additional tags, use 'None' to skip adding the tag.

Version: %s
""" %(VERSION))

parser.add_argument('--inbam', '-i', help= 'Input bam or sam with header, use - to read from stdin (default)', required= True)
parser.add_argument('--aseq', '-aq', help= 'Tag for query alignment sequence. Default "YQ"', default= 'YQ')
parser.add_argument('--alen', '-al', help= 'Tag for reference alignment length. Default "YL"', default= 'YL')
parser.add_argument('--astart', '-as', help= 'Tag for reference start. Default "YS"', default= 'YS')
parser.add_argument('--aend', '-ae', help= 'Tag for reference end. Default "YE"', default= 'YE')
parser.add_argument('--outfmt', '-f', help= 'Output format. Default compressed bam (b)', default= 'b', choices= ['b', 'bu', 'h', 's'])
parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)

args= parser.parse_args()

sam= pysam.AlignmentFile(args.inbam)

if args.outfmt == 's':
    fmt= ''
else:
    fmt= args.outfmt
outheader= sam.header
outheader['PG'].append({'PN': PN, 'ID': PN, 'VN': VERSION, 'CL': ' '.join(sys.argv)})
out= pysam.AlignmentFile('-', 'w%s' %(fmt), header= outheader)

for line in sam:
    if args.aseq != 'None':
        line.set_tag(args.aseq, line.query_alignment_sequence)
    if args.alen != 'None':
        line.set_tag(args.alen, line.reference_length)
    if args.astart != 'None':
        line.set_tag(args.astart, line.reference_start)
    if args.aend != 'None':
        line.set_tag(args.aend, line.reference_end)
    out.write(line)

sam.close()
out.close()
sys.exit()