#!/usr/bin/env python

import sys
import argparse
import re

parser = argparse.ArgumentParser(description= """
DESCRIPTION
    Convert tabular file to table formatted for redmine

EXAMPLE


|_.FILENAME |_. LIBRARY |_. UNPAIRED_READS_EXAMINED |_. READ_PAIRS_EXAMINED |_. UNMAPPED_READS |_. UNPAIRED_READ_DUPLICATES |_. READ_PAIR_DUPLICATES |_. READ_PAIR_OPTICAL_DUPLICATES |_. PERCENT_DUPLICATION |_. ESTIMATED_LIBRARY_SIZE|
|mjb050_E14oxBSAD01 | Unknown Library | 0 | 262372312 | 0 | 0 | 68421371 | 395483 | 0.26078 | 412713349|
|mjb051_E14oxBSAD02 | Unknown Library | 0 | 249082147 | 0 | 0 | 54443524 | 330432 | 0.218577 | 485408946|
|mjb052_E14oxBSAD03 | Unknown Library | 0 | 282625823 | 0 | 0 | 90534858 | 486614 | 0.320335 | 341860496|
|mjb053_E14BSAD04 | Unknown Library | 0 | 237375407 | 0 | 0 | 42147060 | 270606 | 0.177554 | 589634340|
|mjb054_E14BSAD05 | Unknown Library | 0 | 250273348 | 0 | 0 | 49954914 | 306798 | 0.199601 | 542843301|
|mjb055_E14BSAD06 | Unknown Library | 0 | 248324589 | 0 | 0 | 44245629 | 308193 | 0.178177 | 614621734|
   
""", formatter_class= argparse.RawTextHelpFormatter, version= '0.1.0')

parser.add_argument('infile',
                    default= None,
                   help='''Input file to convert. Read from stdin if -
''')

parser.add_argument('--sep', '-s',
                    default= '\s+',
                   help='''Regex to use aS column separator. Default to any blank space ('\\s+'). 
''')

parser.add_argument('--header', '-H',
                    action= 'store_true',
                   help='''First line will be formatted in bold as header.
''')

args= parser.parse_args()

if args.infile is None or args.infile == '-':
    fin= sys.stdin
else:
    fin= open(args.infile, 'rU')

#if args.sep is None:
#    sep= '\t'
#else:
#    sep= args.sep

header= args.header



for line in fin:
    line= re.split(args.sep, line.strip())
    if header:
        line= ' |_. '.join(line)
        line=  '|_. ' + line
        header= False
    else:
        line= '|' + ' | '.join(line)
    print(line + "|")

sys.exit()

