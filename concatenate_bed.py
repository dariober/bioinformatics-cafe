#!/usr/bin/env python

import os
import argparse
import sys
import re
import gzip

parser = argparse.ArgumentParser(description= """

DESCRIPTION:

    Concatenate several bed files in a single one and adds a column identifying
    for each line the source file.
    An error is thrown if one tries to concatenate a file to itself.
    
    This program concatenate any file, not just beds and no checking is done
    about whether the files have bed formatting.
    
EXAMPLES:
    
    Concatenate all bed files in current dir and strip the .bed extension
    
    ls *.bed | concatenate_bed.py -i - -o concoutput.bed -s .bed

    Output (last column is the file name):
--------------------------------------------------------------------------------
chr1    3521705 3521924 CpG:_27 2       52      219     0.2374429       bham359.coverage
chr1    3660700 3661155 CpG:_34 11      302     455     0.6637363       bham359.coverage
chr1    3661735 3662237 CpG:_45 10      311     502     0.6195219       bham359.coverage
...
--------------------------------------------------------------------------------

TODO:
      
""", formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('-i', '--input',
                    nargs='+',
                    type= str,
                    required= True,
                    help="""Input bed files(s) to be concatebated. Use -i - to
read the list of files from stdin (e.g. ls *.bed | concatenate_bed.py -i -)
                    """)

parser.add_argument('-o', '--output',
                    type= str,
                    required= True,
                    help="""Name of output file.
                    """)

parser.add_argument('-s', '--strip',
                    type= str,
                    nargs= '+',
                    required= False,
                    help="""List of regexs to strip from the input file name. Each of these
regexs will be removed. E.g.
-s \.merged\.bed$ \.bed$
will strip both .bed and .merged.bed. Note: Order matters put less specific regex first (that is .merged.bed bofore .bed)
The resulting name will be used as identifier of the source file. Default is not
to strip anything.
                """)

parser.add_argument('-d', '--dir',
                    action= 'store_true',
                    required= False,
                    help="""With this flag the file id in output will include
the directory path. Defualt is to strip the path. 
                """)

parser.add_argument('--skip',
                    type= int,
                    default= 0,
                    required= False,
                    help="""Skip this many lines from each input bed before
writing out (e.g. use --skip 1 to skip the header). 
                """)

parser.add_argument('--fill',
                    type= str,
                    default= None,
                    required= False,
                    help="""If the concatenated bed files have different number of fields,
fill short rows with this string (e.g. --fill NA). Note: The file name will remain
the last column. Default is None meaning don't fill up rows 
                """)

#parser.add_argument('-H', '--header',
#                    type= str,
#                    required= False,
#                    help="""Optional string to use as header.
#                    """)

args = parser.parse_args()

if args.input == ['-']:
    args.input= sys.stdin.readlines()
    args.input= [x.strip() for x in args.input]
if args.output in args.input:
    sys.exit('%s error: Output file %s also present in input list' %(os.path.basename(__file__), args.output))


#if args.header:
#    print(args.header)

fout= open(args.output, 'w')

ncols= 0 ## Keep track of the number of columns in eah line in order to fill in short rows
concbed= []
for f in args.input:
    n= 0
    file_id= f
    if args.dir is False:
        file_id= os.path.split(f)[1]
    if args.strip:
        for r in args.strip:
            file_id= re.sub(r, '', file_id)
    print('Concatenating: %s; File ID: %s' %(f, file_id))

    if f.endswith('.gz')
        fin= gzip.open(f)
    else:
        fin= open(f)
        
    for line in fin:
        if n < args.skip:
            n += 1
            continue
        line= line.rstrip('\n\r')
        line= line.split('\t')
        if len(line) > ncols:
            ncols= len(line)
        line.append(file_id)
        concbed.append(line)
    fin.close()

ncols= ncols + 1 ## +1 because file name has been appended to each row

if args.fill is not None:
    " Memo: concbed is a list of list. Each inner list a bed row "
    for i in range(0, len(concbed)):
        line= concbed[i]
        if len(line) < ncols: 
            fill= [args.fill] * (ncols - len(line))
            line= line[:-1] + fill + [line[-1]]
            concbed[i]= line

for line in concbed:
    fout.write('\t'.join(line) + '\n')
fout.close()
