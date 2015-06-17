#!/usr/bin/env python

import argparse
import sys
import glob
import os
import gzip
import re

parser = argparse.ArgumentParser(description= """
DESCRIPTION
    Concatenate files and add a column of file identifiers. The output is
    going to be in "database long format".
    
    Typical use case: You have a pipeline that generates file in tabular format,
    one file per input, all files the same format. You want to concatenate these
    tables while keeping the identity of the files and optionally add une header
    for the cat'd files.

EXAMPLE
    cd tableCat
    ./tableCat.py -i test/*.txt -H
""", formatter_class= argparse.RawTextHelpFormatter, prog= os.path.basename(__file__))

parser.add_argument('--input', '-i',
    required= True,
    nargs= '+',
    help='''Input files to concatenate
''')

parser.add_argument('--sep', '-s',
    default= '\t',
    help='''Separator used to right-bind the column with file id. Default to tab.
''')

parser.add_argument('--keepdir', '-kd',
    action= 'store_true',
    help='''Keep the the directory name in the output. Default is to strip it.
''')

parser.add_argument('--regex', '-r',
    default= None,
    help='''an optional regex to strip from the file name to make the ID nicer.
''')

parser.add_argument('--idColName', '-id',
    default= 'file_id',
    help='''Name for the additional column with the file id. Default 'file_id'.
''')

parser.add_argument('--skip', '-S',
    type= int,
    default= 0,
    help='''Skip this many lines from *each* file.
''')

parser.add_argument('--firsthdr', '-H',
    action= 'store_true',
    help='''First non-skipped line from *first* file is header and it is printed.
    The first non-skipped line from following files is NOT printed.
''')


parser.add_argument('--version', action='version', version='%(prog)s 0.1.0')

# -----------------------------------------------------------------------------

def globToList(glist):
    """Expand glob(s) to list of individual files. Duplicates removed!
    """
    if type(glist) != list:
        sys.exit('Invalid arg')
    globbed= []
    for g in glist:
        xglobbed= glob.glob(g)
        for x in xglobbed:
            globbed.append(x)
    globbed= sorted(list(set(globbed)))
    return globbed

def xopen(x):
    """Seemlessly open flat or gzip file.
    x: File name
    Return: Open file handle ready to iterate through
    """
    try:
        fin= gzip.open(x, 'rb')
        fin.readline()
        fin.seek(0)
    except IOError:
        try:
            fin= open(x)
        except IOError:
            sys.exit('Cannot open %s' %(x))
    return fin

def parseID(filename, keepdir= False, regex= None):
    """Parse the filename to optionally remove dir and a regex
    """
    id= filename
    if not keepdir:
        id= os.path.basename(filename)
    if regex:
        id= re.sub(regex, '', id)
    return(id)

if __name__ == '__main__':
    args= parser.parse_args()
    files= globToList(args.input)
    if files == []:
        sys.exit('No file found to concatenate!')
    is_first_file= True
    for f in files:
        # sys.stderr.write('File: %s ...' %(f))
        fin= xopen(f)
        id= parseID(f, args.keepdir, args.regex)
        is_first_line= True
        n= 0
        skip= args.skip
        for line in fin:
            line2= line.rstrip('\n') + args.sep + id
            if skip > 0:
                skip-=1
                continue
            if args.firsthdr:
                if is_first_file and is_first_line:
                    print line.rstrip('\n') + args.sep + args.idColName
                    n+=1
                else:
                    if is_first_line:
                        pass
                    else:
                        print line2
                        n+=1
                is_first_line= False
            else:
                print line2
                n+=1
        is_first_file= False
        # sys.stderr.write(' %s lines printed\n' %(n))
    sys.exit()
    
