#!/usr/bin/env python3

import argparse
import sys
import subprocess as sp
import os
import distutils.spawn

VERSION= '0.1.0'
thisprog= '%s %s' %(os.path.basename(__file__), VERSION)

parser = argparse.ArgumentParser(description= """
DESCRIPTION
Similar to Unix `split`, split a file in multiple files each sub-file
containing `-l` number of lines at most. Output files are compressed, this the
main difference with Unix split.
""", formatter_class= argparse.RawTextHelpFormatter, prog= os.path.basename(__file__))

parser.add_argument('--infile', '-i',
                   help= '''Input file to split or - to read from stdin. Default %(default)s
                   ''',
                   default='-')

parser.add_argument('--lines', '-l',
                   help= '''Put this many lines per output file. If <= 0, do not split at all. Default %(default)s          
                   ''',
                   type= int,
                   default= -1)

parser.add_argument('--prefix', '-p',
                   help= '''Prefix for output files. Default %(default)s          
                   ''',
                   default= 'x')

parser.add_argument('--suffix_length', '-a',
                   help= '''Suffix length. Suffix characters are numeric. Default %(default)s          
                   ''',
                   type= int,
                   default= 4)

parser.add_argument('--extension', '-x',
                   help= '''Extension to add to output files. Default %(default)s          
                   ''',
                   default= '.gz')

parser.add_argument('--version', '-v', action='version', version='%(prog)s ' + VERSION)

args= parser.parse_args()

if distutils.spawn.find_executable("pigz") != '':
    xzip= 'pigz'
elif distutils.spawn.find_executable("gzip") != '':
    xzip='gzip'
else:
    sys.stderr.write('\nCannot find a suitable system command to compress\n')
    sys.exit(1)

if args.infile == '-':
    fin= sys.stdin
else:
    fin= open(args.infile)

nlines= 0
nprefix= 0
suffix= str(nprefix).zfill(args.suffix_length)
cmd= '%s > %s%s%s' %(xzip, args.prefix, suffix, args.extension)
p= sp.Popen(cmd, stdout= sp.PIPE, stdin= sp.PIPE, universal_newlines= True, shell= True)
for line in fin:
    if args.lines <= 0 or nlines < args.lines:
        p.stdin.write(line)
    else:
        p.stdin.close()
        p.wait()
        nlines= 0
        nprefix += 1
        suffix= str(nprefix).zfill(args.suffix_length)
        cmd= '%s > %s%s%s' %(xzip, args.prefix, suffix, args.extension)
        p= sp.Popen(cmd, stdout= sp.PIPE, stdin= sp.PIPE, universal_newlines= True, shell= True)
        p.stdin.write(line)
    nlines += 1
fin.close()
sys.exit()
