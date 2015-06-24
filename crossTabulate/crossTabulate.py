#!/usr/bin/env python

import os
import sys
import argparse
import re
import gzip

VERSION= '0.2.0'

parser = argparse.ArgumentParser(description= """
DESCRIPTION
    Cross tabulate files. Typical use case: You have a number of tabular files, 
    like bed files, with a column of features and a column of values. You want 
    to make a table where first column is feature name and following columns
    are the value columns from all files.   

""", formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('infiles',
                   nargs= '+',
                   help='''One or more input files to cross-tabulate. Use - to
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
parser.add_argument('--sep', '-sep', type= str, default= '\t', help= 'Field separator for input and output files. Default to tab')

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

if args.infiles == ['-']:
    args.infiles= sys.stdin.readlines()
    args.infiles= [x.strip() for x in args.infiles]
    
if args.keepdir is False:
    header= [os.path.split(x)[1] for x in args.infiles] ## Remove path
else:
    header= [x for x in args.infiles]
if args.strip is not None:
    header= [re.sub(args.strip, '', x) for x in header]
header= ['feature_id'] + header

print(args.sep.join(header))
    
infiles= []
for lib in args.infiles:
    if(lib.endswith('.gz')):
        infiles.append(gzip.open(lib, 'rb'))
    else:
        infiles.append(open(lib))

while True:
    for i in range(0, len(infiles)):
        if i == 0:
            # First file: include gene name and count
            fline= infiles[i].readline().strip().split(args.sep)
            if fline == ['']:
                break
            feature_id= fline[name_col]        
            line= [fline[name_col], fline[cnt_col]]
        else:
            # Following files: Include only count
            fline= infiles[i].readline().strip().split(args.sep)
            if fline == ['']:
                break
            this_feature= fline[name_col]
            if this_feature != feature_id:
                sys.exit('Error: Input files have the same features and/or they are not in the same order: %s, %s' %(feature_id, this_feature))
            line.append(fline[cnt_col])
    if fline == ['']:
        break
    print(args.sep.join(line))
for f in infiles:
    f.close()    
sys.exit()
