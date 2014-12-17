#!/usr/bin/env python

import sys
import pysam

docstring= """
DESCRIPTION
Add to the BAM header the RG tags found in the reads and write out a copy of
the input file with the header updated.

Typical use case: You run `samtools merge -r out.bam in1.bam in2.bam` and you
get out.bam with reads tagged with RG. However, out.bam does not have those
RG tags in the header. Fix out.bam with this script.

USAGE
addRGtoSAMHeader.py <in.bam> <out.bam>

Use - to write BAM to stdout.

IMPORTANT:
* Additional attributes in @RG (SM, PL, ...) are not included.
* If @RG lines are found in the header, the new RG tags are appended.
* Read w/o RG tag go to output as they are, w/o warnings.
"""

if len(sys.argv) != 3 or sys.argv[1] in ('-h', '--help'):
    print docstring
    sys.exit(1)

# ------------------------------------------------------------------------------
def rgSetToListOfDict(RGset):
    """Convert the set of RG tags to a list of dictionaries ready to be add
    to the original bam header.
    """
    dlist= []
    for rg in RGset:
        dlist.append({'ID': rg})
    return dlist

def addRGsetToHeader(header, RGset):
    """Return the header dictionary with the additional RG entries included.
    If RG tags are found in header, the new ones are appended.
    
    header:
        comes from the input bam file
    RGseq:
        Typically from rgSetToListOfDict()

    Output is dictionary header ready for pysam.AlignmentFile(..., header= header)
    """
    header= insam.header
    outheader= {}
    for k in header.keys():
        outheader[k]= header[k]
    
    rg= rgSetToListOfDict(RGset)
    if 'RG' in outheader:
        current= outheader['RG']
        for x in rg:
            if x not in current:
                current.append(x)
        outheader['RG']= current
    else:
        outheader['RG']= rg
    return outheader

# ------------------------------------------------------------------------------

infile= sys.argv[1]
outfile= sys.argv[2]

# ------------------------------------------------------------------------------

insam= pysam.AlignmentFile(infile, "rb")

## Collect RG tags
sys.stderr.write('Collecting RG tags from %s... ' %(infile))
i= 0;
RGset= set()
for aln in insam:
    try:
        RGset.add(aln.opt('RG'))
    except KeyError:
        pass
    i += 1
sys.stderr.write('%s tags; %s reads\n' %(len(RGset), i))


## Prepare output file:
outheader= addRGsetToHeader(insam.header, RGset)
outsam= pysam.AlignmentFile(outfile, "wb", header= outheader)

## Write out to new file:
insam.close()
insam= pysam.AlignmentFile(infile, "rb")
sys.stderr.write('Copying reads to %s\n' %(outfile))
j= 0
for aln in insam:
    j+=1
    outsam.write(aln)
insam.close()
outsam.close()

## Minimal check:
if(i != j):
    sys.stderr.write('\nERROR: Reads in input (%s) != read in output (%s)!\n\n' %(i, j))
    sys.exit(1)
    
sys.exit()
