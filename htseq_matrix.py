#!/home/berald01/.local/bin/python

import os
import sys
import argparse
import re

VERSION= '0.1.0'

parser = argparse.ArgumentParser(description= """
DESCRIPTION
    Produces a matrix of genes (rows) and libraries (columns) from the output of
    htseq-count. Output will be suitable for edgeR/DESeq after having removed 
    unwanted rows like "__no_feature"
    

EXAMPLE
    # Matrix of all files ending in .htseq, remove suffix. 
    # Exclude rows beginning with '__'. 
    htseq_matrix.py *.htseq -s '\.htseq$' | grep -v '^__'

""", formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('htseqfiles',
                   nargs= '+',
                   help='''One or more output files from htseq-count to be converted
to matrix. Output files should have the same features in the same order. Use - to
read the list of files from stdin.
                   ''')

parser.add_argument('--strip', '-s',
                   required= False,
                   type= str,
                   default= None,
                   help='''Remove from each file name this regular expression.
The replacement is applied to the file name only unless --keepdir is active. With
--keepdir the replacement is applied to the whole path as passed to the list of files.
                   ''')

parser.add_argument('--keepdir', '-d',
                   required= False,
                   action= 'store_true',
                   help='''Do not remove directory path from file names (by default
it is removed).
                   ''')
parser.add_argument('--version', '-v', action='version', version='%(prog)s %(VERSION)s' %{'prog': '%(prog)s', 'VERSION': VERSION})

args = parser.parse_args()

# ------------------------------------------------------------------------------

NAME_COL= 0 # Index of feature name and count columns
CNT_COL= 1

if args.htseqfiles == ['-']:
    args.htseqfiles= sys.stdin.readlines()
    args.htseqfiles= [x.strip() for x in args.htseqfiles]
    
if args.keepdir is False:
    header= [os.path.split(x)[1] for x in args.htseqfiles] ## Remove path
else:
    header= [x for x in args.htseqfiles]
if args.strip is not None:
    header= [re.sub(args.strip, '', x) for x in header]
header= ['feature_id'] + header

print('\t'.join(header))
    
infiles= []
for lib in args.htseqfiles:
    infiles.append(open(lib))

while True:
    for i in range(0, len(infiles)):
        if i == 0:
            # First file: include gene name and count
            fline= infiles[i].readline().strip().split('\t')
            feature_id= fline[NAME_COL]
            if fline == ['']:
                break
            line= [fline[NAME_COL], fline[CNT_COL]]
        else:
            # Following files: Include only count
            fline= infiles[i].readline().strip().split('\t')
            if fline == ['']:
                break
            this_feature= fline[NAME_COL]
            if this_feature != feature_id:
                sys.exit('Input files do not seem to have the same features and/or they are not in the same order: %s, %s' %(feature_id, this_feature))
            line.append(fline[CNT_COL])
    if fline == ['']:
        break
    print('\t'.join(line))
for f in infiles:
    f.close()    
sys.exit()
