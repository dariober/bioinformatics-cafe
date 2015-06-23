#!/usr/bin/env python

import os
import sys
import argparse
import re
import gzip

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
read the list of files from stdin. File ending with .gz are opened with gzip.
                   ''')

parser.add_argument('--strip', '-s',
                   required= False,
                   type= str,
                   default= None,
                   help='''Remove from each file name this regular expression.
The replacement is applied to the file name only unless --keepdir is active. With
--keepdir the replacement is applied to the whole path as passed to the list of files.
                   ''')

parser.add_argument('--nameCol', '-n', type= int, default= 1, help= 'Column index of the gene or feature name. 1-based, default 1 (1st column)')
parser.add_argument('--cntCol', '-c', type= int, default= 2, help= 'Column index of the counts. 1-based, default 2 (2nd column)')

parser.add_argument('--keepdir', '-d',
                   required= False,
                   action= 'store_true',
                   help='''Do not remove directory path from file names (by default
it is removed).
                   ''')

parser.add_argument('--version', '-v', action='version', version='%(prog)s %(VERSION)s' %{'prog': '%(prog)s', 'VERSION': VERSION})

args = parser.parse_args()

# ------------------------------------------------------------------------------

name_col= args.nameCol - 1
cnt_col= args.cntCol - 1

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
    if(lib.endswith('.gz')):
	infiles.append(gzip.open(lib, 'rb'))
    else:
        infiles.append(open(lib))

while True:
    for i in range(0, len(infiles)):
        if i == 0:
            # First file: include gene name and count
            fline= infiles[i].readline().strip().split('\t')
            if fline == ['']:
                break
            feature_id= fline[name_col]        
            line= [fline[name_col], fline[cnt_col]]
        else:
            # Following files: Include only count
            fline= infiles[i].readline().strip().split('\t')
            if fline == ['']:
                break
            this_feature= fline[name_col]
            if this_feature != feature_id:
                sys.exit('Input files do not seem to have the same features and/or they are not in the same order: %s, %s' %(feature_id, this_feature))
            line.append(fline[cnt_col])
    if fline == ['']:
        break
    print('\t'.join(line))
for f in infiles:
    f.close()    
sys.exit()
