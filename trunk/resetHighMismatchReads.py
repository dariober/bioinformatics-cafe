#!/usr/bin/env python

import sys
import pysam
import argparse
parser = argparse.ArgumentParser(description= """
DESCRIPTION
    Set mapping quality to 0 for reads with high density of mismatches and reset
    read quality to "!".
 
""", formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('--input', '-i',
                   required= True,
                   help='''Input sam or bam file or - to read *bam* from stdin

''')

parser.add_argument('--output', '-o',
                   required= False,
                   default= '-',
                   help='''Output sam or bam file, format determined by extension. Default is sam to
stdout. SAM includes header.

''')

parser.add_argument('--mismDensity', '-m',
                   required= False,
                   type= float,
                   default= 0.1,
                   help='''Max allowed density of mismatches expressed as NM / aligned length. E.g
0.1 allows 1 mismatch in a 10bp alignment. Default 0.1.
''')

parser.add_argument('--mtag', '-mt',
                   required= False,
                   default= 'YM',
                   help='''Original mapq score stored will be stored in this tag.
Existing tag will be overwritten. Default YM.
''')

parser.add_argument('--qtag', '-qt',
                   required= False,
                   default= 'YQ',
                   help='''Original quality string will be stored in this tag.
Existing tag will be overwritten. Default YQ.
''')


parser.add_argument('--tag', '-t',
                   required= False,
                   default= 'NM',
                   help='''Tag to use to extract the number of mismatches. Default 'NM'
Tag from BSSeqMismatches.jar is XZ. If the tag is not found the record is written
unchanged.
''')

args= parser.parse_args()
# -----------------------------------------------------------------------------

def getTagValue(tags, tag):
    """Returns the value associated to tag given the list of tuples tags
    """
    tv= [x[1] for x in tags if x[0] == tag]    
    if (len(tv) > 1):
        sys.stderr.write("More than one value found. Returning first one\n")
    if tv == []:
        return(None)
    else:
        return(tv[0])
    
# ------------------------------------------------------------------------------

# Input
samfile = pysam.Samfile(args.input)

if args.output.endswith('.sam') or args.output == "-":
    outsam= pysam.Samfile(args.output, "wh", template= samfile)
elif args.output.endswith('.bam'):
    outsam= pysam.Samfile(args.output, "wb", template= samfile)
else:
    sys.exit("Invalid extension for output file. Must be .sam or .bam. Got: " + args.output)

## N. mismatches
mm= args.mismDensity
if mm < 0 or mm > 1:
    sys.exit("Invalid density of mismatches. Must be between 0 and 1, got " + str(mm))

# -----------------------------------------------------------------------------

for aln in samfile:
    NM= getTagValue(aln.tags, args.tag)
    if (aln.alen is not None) and (aln.alen > 0) and (NM is not None):
        if (NM / float(aln.alen)) > mm:   
            aln.tags= aln.tags + [(args.mtag, int(aln.mapq))] ## Put MAPQ in its tag value and reset
            aln.mapq= 0
            aln.tags= aln.tags + [(args.qtag, aln.qual)] ## Put quality string  in its tag and reset
            aln.qual= "!" * len(aln.qual)
    outsam.write(aln)

samfile.close()
outsam.close()