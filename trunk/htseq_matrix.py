#!/home/berald01/.local/bin/python

import os
import sys
import argparse
import re
import redmine

parser = argparse.ArgumentParser(description= """
DESCRIPTION
    Produces a matrix of genes (rows) and libraries (columns) from the output of
    htseq-count. Output will be suitable for edgeR/DESeq.
    Stdout will print the rows from "no_feature" to the end.

EXAMPLE
    ## Produce matrix of all the files ending in .htseq, remove .htseq from file name
    find . -name '*.htseq' | sort | htseq_matrix.py - -s '\.htseq$' htseq.matrix
TODO:

""", formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('htseqfiles',
                   nargs= '+',
                   help='''One or more output files from htseq-count to be converted
to matrix. Output files should have the same features in the same order. Use - to
read the list of files from stdin.
                   ''')

parser.add_argument('--outmatrix', '-o',
                   required= True,
                   type= str,
                   help='''Output file where the matrix will be written to.
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

parser.add_argument('--redmine',
                    action= 'store_true',
                   help='''Format the report output as a table for redmine (default
is a tab separatetd table).
                   ''')

args = parser.parse_args()

# ------------------------------------------------------------------------------

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

outf= open(args.outmatrix, 'w')
outf.write('\t'.join(header) + '\n')
if args.redmine:
    print(redmine.list2table(header, header= True))
else:
    header= '\t'.join(header)
    print(header)
    
is_feature= True ## Mark to say whether the current line is feature to be inlcuded in the matrix
                 ## Will turn to False once the row 'no_feature' is found

infiles= []
for lib in args.htseqfiles:
    infiles.append(open(lib))
read_count= [0] * len(infiles) ## Keep track of how many reads are mapped to feature

while True:
    for i in range(0, len(infiles)):
        if i == 0:
            fline= infiles[i].readline()
            fline= fline.strip().split('\t')
            feature_id= fline[0]
            if feature_id == 'no_feature':
                is_feature= False
            if fline == ['']:
                break
            line= fline
            if is_feature:
                read_count[i] += int(line[1])
        else:
            fline= infiles[i].readline()
            fline= fline.strip().split('\t')
            if fline == ['']:
                break
            next_feature= fline[0]
            if next_feature != feature_id:
                sys.exit('Input files do not seem to have the same features and/or they are not in the same order: %s, %s' %(feature_id, next_feature))
            line.append(fline[1])
            if is_feature:
                read_count[i] += int(fline[1])
    if fline == ['']:
        break
    if is_feature:
        outf.write('\t'.join(line) + '\n')
    else:
        if args.redmine:
            print(redmine.list2table(line))
        else:
            print('\t'.join(line))
read_count= [str(x) for x in read_count]
if args.redmine:
    print(redmine.list2table(['assigned'] + read_count))
else:
    print('\t'.join(['assigned'] + read_count))
for f in infiles:
    f.close()    
outf.close()
sys.exit()