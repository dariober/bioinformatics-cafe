#!/usr/bin/env python

import sys
import pysam
import argparse
import os
# import pprint <- for debugging useful to print list and dict in readable form

VERSION= '0.1.0'
thisprog= '%s %s' %(os.path.basename(__file__), VERSION)

parser = argparse.ArgumentParser(description= """
DESCRIPTION
Add to the BAM header the RG tags found in the reads and write out a copy of
the input file with the header updated.

Typical use case: You run `samtools merge -r out.bam in1.bam in2.bam` and you
get out.bam with reads tagged with RG. However, out.bam does not have those
RG tags in the header. Fix out.bam with this script.

IMPORTANT:
* A pre-existing @RG key in input is completely replaced by the new one!
* Read w/o RG tag go to output as they are, w/o warnings.

Requires: pysam 0.8.1
""", formatter_class= argparse.RawTextHelpFormatter, prog= os.path.basename(__file__))

parser.add_argument('--input', '-i',
                   required= True,
                   help='''Input bam file.
                   ''')

parser.add_argument('--out', '-o',
                   required= True,
                   help='''Output bam file. Use - to write to stdout
                   ''')

parser.add_argument('--rglb', '-lb',
                   default= 'NA',
                   help= '''Read Group Library. Default: NA
                   ''')

parser.add_argument('--rgsm', '-sm',
                   default= 'NA',
                   help= '''Read Group sample name. Default: NA
                   ''')

parser.add_argument('--rgpl', '-pl',
                   default= 'NA',
                   help= '''Read Group platform (e.g. illumina, solid). Default: NA
                   ''')

parser.add_argument('--rgpu', '-pu',
                   default= 'NA',
                   help= '''Read Group platform unit (eg. run barcode). Default: NA
                   ''')

parser.add_argument('--rgpg', '-pg',
                   required= False,
                   help= '''Programs used for processing the read group. Default "%s"
                   ''' %(thisprog))

parser.add_argument('--onlyhdr', '-H',
                   action= 'store_true',
                   help='''Output only the modified header, w/o appending reads.
                   ''')

parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)

args= parser.parse_args()

"""
RGCN (String)	Read Group sequencing center name Default value: null.
RGDS (String)	Read Group description Default value: null.
RGDT (Iso8601Date)	Read Group run date Default value: null.
RGPI (Integer)	Read Group predicted insert size Default value: null. 
"""

# ------------------------------------------------------------------------------
def rgSetToListOfDict(RGset, LB, SM, PL, PU, PG):
    """Convert the set of RG tags to a list of dictionaries ready to be add
    to the original bam header.
    
    Output list will look like:
    [{'ID': 'myid1',
         'SM': 'test' ...},
        {'ID': 'myid2',
         'LB': 'test',
         'SM': 'test' ...},
        {'ID': 'myid3',
         'SM': 'test' ...}]    
    """
    dlist= []
    for rg in RGset:
        dlist.append({'ID': rg, 'LB': LB, 'SM': SM, 'PL': PL, 'PU': PU, 'PG': PG})
    return dlist

def addRGsetToHeader(header, RGlist):
    """Return the header dictionary with the additional RG entries included.
    If RG tags are found in header, the new ones are appended.
    
    header:
        comes from the input bam file
    RGseq:
        List of dictionary. Typically from rgSetToListOfDict()

    Output is dictionary header ready for
    pysam.AlignmentFile(..., header= header)
    
    MEMO: The sam header in pysam is dict looking like this:
    {'SQ': [{'LN': 16571, 'SN': 'chrM'},
            {'LN': 249250621, 'SN': 'chr1'},
            {'LN': 59373566, 'SN': 'chrY'}],
     'PG': [{'PN': 'bwa', 'ID': 'bwa', 'VN': '0.7.10-r789'},
            {'PN': 'MarkDuplicates', 'ID': 'MarkDuplicates', 'VN': '1.127'}]}
    """
    header= insam.header
    outheader= {}
    for k in header.keys():
        outheader[k]= header[k]
    outheader['RG']= RGlist
    return outheader

# ------------------------------------------------------------------------------

insam= pysam.AlignmentFile(args.input, "rb")

## Collect RG tags
sys.stderr.write('Collecting RG tags from %s... ' %(args.input))
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
if not args.rgpg:
    rgpg= thisprog
else:
    rgpg= args.rgpg
    
RGlist= rgSetToListOfDict(RGset, LB= args.rglb, SM= args.rgsm, PL= args.rgpl, PU= args.rgpu, PG= rgpg)

outheader= addRGsetToHeader(insam.header, RGlist)
outsam= pysam.AlignmentFile(args.out, "wb", header= outheader)

if args.onlyhdr:
    outsam.close()
    sys.exit()

## Write out to new file:
insam.close()
insam= pysam.AlignmentFile(args.input, "rb")
sys.stderr.write('Copying reads to %s\n' %(args.out))
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
