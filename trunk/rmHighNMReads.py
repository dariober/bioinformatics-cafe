#!/usr/bin/env python

import sys
import pysam

docstring= """DESCRIPTION

Set mapping quality to 0 for reads with high density of mismatches.

USAGE

rmHighNMreads.py <in.bam> <density mism> [out.bam|.sam]

<in.bam>:
    Input sam or bam file or - to read *bam* from stdin

<density mism>:
    Max allowed density of mismatches expressed as NM / aligned length. E.g
    0.1 allows 1 mismatch in a 10bp alignment.

[out.bam|.sam]:
    Output sam or bam file, format determined by extension. Defualt is sam to
    stdout. SAM includes header.
"""

def getTagValue(tags, tag):
    """Returns the value associated to tag given the list of tuples tags
    """
    tv= [x[1] for x in tags if x[0] == tag]
    if (len(tv) > 1):
        sys.stderr.write("More than one value found. Returning first one\n")
    return(tv[0])
    
# ------------------------------------------------------------------------------
if (len(sys.argv) not in [3, 4] or sys.argv[1] in ['-h', '--help']):
    print(docstring)
    sys.exit(1)

# Parse args -------------------------------------------------------------------
# Input
inbam= sys.argv[1]
samfile = pysam.Samfile(inbam)

# Output
if len(sys.argv) == 4:
    outfile= sys.argv[4]
    if outfile.endswith('.sam'):
        outsam= pysam.Samfile(outfile, "wh", template= samfile)
    elif outfile.endswith('.bam'):
        outsam= pysam.Samfile(outfile, "wb", template= samfile)
    else:
        sys.exit("Invalid extension for output file. Must be .sam or .bam. Got: " + outfile)
else:
    outsam= pysam.Samfile("-", "wh", template= samfile)

## N. mismatches
mm= float(sys.argv[2])
if mm < 0 or mm > 1:
    sys.exit("Invalid density of mismatches. Must be between 0 and 1, got " + str(mm))

# -----------------------------------------------------------------------------

for aln in samfile:
    NM= getTagValue(aln.tags, 'NM')
    L= float(aln.alen)
    if (L > 0) and (NM / L > mm):
        aln.mapq= 0
    outsam.write(aln)

samfile.close()
outsam.close()