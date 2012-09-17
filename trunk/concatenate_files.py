#!/home/berald01/.local/bin/python

import sys
import argparse
import re

parser = argparse.ArgumentParser(description= """
DESCRIPTION
    Concatenate files and prepened to each row the file name (tab separated), optionally strip
    a pyhton regexp from the file name.

       
EXAMPLE
    ## Concatenate all files ending in .txt, strip the .txt extension and skip first line from each file:
    ls *.txt | concatenate_files.py --skip 1 --strip '\.txt$' - > mycat.txt
DEPENDS-ON:
    
TODO:

""", formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('filelist',
                   nargs= '+',
                   help='''List of files to concatenate, use - to read this list from stdin
                   ''')

parser.add_argument('-s', '--strip',
                   type= str,
                   default= None,
                   required= False,
                   help='''Strip this regular expression from each filename.
                   ''')

parser.add_argument('-S', '--skip',
                   type= int,
                   default= 0,
                   required= False,
                   help='''Skip this many lines from each input file (default 0)
                   ''')

args = parser.parse_args()

# -----------------------------------------------------------------------------

if args.filelist == ['-']:
    args.filelist= sys.stdin.readlines()
    args.filelist= [x.strip() for x in args.filelist]

for f in args.filelist:
    n= 0
    fin= open(f)
    fname= f
    if args.strip is not None:
        fname= re.sub(args.strip, '', fname)
    while n < args.skip:
        fin.readline()
        n += 1
    for line in fin:
        print(fname + '\t' + line.rstrip('\r\n'))
    fin.close()

sys.exit()