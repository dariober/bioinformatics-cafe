#!/usr/bin/env python

import sys
import argparse
import pysam

parser = argparse.ArgumentParser(description= """
DESCRIPTION
    Reformat the SAM output from BSmooth/bswc_bowtie2_align.pl to replace the
    sequence in the SEQ field with the sequene in the sequence in YO tag.
    
    The YO tag is reset to '' and a new tag, YB, will contain the converted sequence
    (now replaced by the sequence in YO)
    
    This is done because (from bswc_bowtie2_align.pl):
        YO:Z:   The original read sequence.  The SAM SEQ field contains the
        C-to-T-converted sequence.
    Conversion is required for methylation calling with mpileup2methylation.py

USAGE
    swapBSmoothSAM.py -i <myreads.sam> > swapped.sam

EXAMPLE:
    ## Read from bam, output bam
    samtools view -b to017_7_bs_epi.bam | swapBSoothSAM.py -b -i -
    
    ## Read and write sam Note the -h option to samtools
    samtools view -h to017_7_bs_epi.bam | swapBSoothSAM.py -S -i -

MEMO: For Defaults and flags are consistent with samtools:
    Default input: bam (only relevant for stdin),
    Default output: sam
    -b: write bam
    -S: input is sam (only relevant for stdin)
""", formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('--input', '-i',
                   required= True,
                   help='''Input sam file or - to read from stdin
                   ''')

parser.add_argument('--isam', '-S',
                   required= False,
                   action= 'store_true',
                   help='''Input stream is SAM. Default is BAM.
                   ''')

parser.add_argument('--outbam', '-b',
                   required= False,
                   action= 'store_true',
                   help='''Write to stdout as BAM. Defualt is to write SAM.
                   ''')

args = parser.parse_args()

# -----------------------------------------------------------------------------

def getTag(tagList, tag):
    """Return the index associated to to the tag from the list
    of tags in pysam.AlignedRead.tags
    E.g. getTag([('AS', -8), ('XN', 0), ('XM', 0)], 'AS') >>> 0
    """
    i= 0
    for t, v in tagList:
        if t == tag:
            return(i)
        i += 1
    return(None)

def revcomp(x):
    """Reverse complement DNA sequence x
    x= 'ACTGN'
    revcomp(x) 
    """
    revDict= {'A': 'T',
              'C': 'G',
              'G': 'C',
              'T': 'A',
              'N': 'N'}
    xrc= []
    xrc= [revDict[n] for n in x]
    xrc= xrc[::-1]
    return(''.join(xrc))
   
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
    yo_index= getTag(read.tags, 'YO')
    readSeq= read.tags[yo_index][1]
    if readSeq is None:
        sys.exit('YO tag not found on line %s' %(read))
    strand= read.tags[getTag(read.tags, 'XB')][1] ## Get strand: Watson (+) or Crick (-)
    if strand == 'C':
        readSeq= revcomp(readSeq)
        read.flag += 16
    elif strand == 'W':
        pass
    elif strand in ['WR', 'CR']:
        sys.exit('Alignment to WR or CR not supported!')
    else:
        sys.exit('Invalid XB tag or XB tag not found')
    convSeq= read.seq
    read.seq= readSeq
    newTagList= read.tags
    newTagList[yo_index]= ('YO', '')
    newTagList= newTagList + [("YB", convSeq)]
    read.tags= newTagList
    outfile.write(read)
samfile.close()
outfile.close()

sys.exit()
