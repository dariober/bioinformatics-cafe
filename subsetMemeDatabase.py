#!/usr/bin/env python

import argparse
import sys

parser= argparse.ArgumentParser()

argparse.ArgumentParser(description= """
DESCRIPTION:
    Select motifs from meme database to create a database subset.

Memo: The line to parse to match motifs names looks like:
MOTIF MA0004.1 Arnt

EXAMPLE
    cat mot.txt
    MA0113.2
    MA0012.1
    MA0015.1
    MA0004.1

    subsetMemeDatabase.py -m mot.txt -db JASPAR_CORE_2014.meme

    ## Extract by TFBS common name.
    ## E.g. all the motifs matching FKH
    grep -P -i 'MOTIF .*fkh.*' JASPAR_CORE_2014.meme \
    | cut -d ' ' -f 2 \
    | subsetMemeDatabase.py -m - -db JASPAR_CORE_2014.meme

""", formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('--memedb', '-db',
                   required= True,
                   help='''Input meme database to parse. 
                   ''')

parser.add_argument('--motifs', '-m',
                   required= True,
                   help='''File with meme motif IDs to extract. One per line.
The motif ID looks like 'MA0004.1'. Use - to read from stdin. 
                   ''')

args= parser.parse_args()

# ------------------------------------------------------------------------------

if args.motifs == '-':
    inmot= sys.stdin
else:
    inmot= open(args.motifs)

motifs= []
for line in inmot:
    motifs.append(line.strip())

L= len(motifs)
motifs= list(set(motifs))
if L != len(motifs):
    sys.stderr.write('Warning: Duplicate motif IDs found in input motifs\n')

inmot.close()

indb= open(args.memedb)

header= True
for line in indb:
    line= line.strip()
    ## Print lines as they appear in meme db untile the first MOTIF is found
    if line.startswith('MOTIF'):
        motif= line.split()[1]
        header= False
    if header:
        print(line)
        continue
       
    if motif in motifs:
        print(line)
