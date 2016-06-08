#!/usr/bin/env python

import sys
import argparse
import gzip

VERSION= '0.2.0'

parser = argparse.ArgumentParser(description= """
DESCRIPTION
    Extract DNA sequences into a fasta file based on feature coordinates.
    
    Same jobs as bedtools fastaFromBed with the difference that the fasta file is loaded
    in memory. This makes fastaFromBed.py much faster than bedtools when extracting
    several sequences (e.g. millions).
    
    See also 
   
""", formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('--fasta', '-fi',
                   help='''Input fasta file. Can be gzipped (but slower)

                   ''')

parser.add_argument('--bed', '-bed',
                   help='''BED file of ranges to extract from -fi. Use - to read
from stdin.

                   ''')

parser.add_argument('--fo', '-fo',
                   default= None,
                   required= True,
                   help='''Output file (can be FASTA or TAB-delimited). Use -
for stdout.
                    ''')

parser.add_argument('--name', '-name',
                    action= 'store_true',
                   help='''Use the name field for the FASTA header.
        
                   ''')

parser.add_argument('--tab', '-tab',
                    action= 'store_true',
                   help='''Write output in TAB delimited format. Default is FASTA format.
       
                   ''')

parser.add_argument('--strand', '-s',
                    action= 'store_true',
                   help='''Force strandedness. If the feature occupies the antisense,
strand, the sequence will be reverse complemented.
- By default, strand information is ignored.
                   ''')

parser.add_argument('--verbose', '-V',
                    action= 'store_true',
                   help='''Show some progress report
                   ''')

parser.add_argument('--version', action='version', version='%s' %(VERSION))
   
args= parser.parse_args()

def DNAreverseComplement(x):
    """Reverse complemnt DNA string.
    Converison table from http://arep.med.harvard.edu/labgc/adnan/projects/Utilities/revcomp.html
    DNAreverseComplement('ACTGNactgn') -> ncagtNCAGT
    """
    rcDict= {
        'A': 'T',
        'T': 'A',
        'U': 'A',
        'G': 'C',
        'C': 'G',
        'Y': 'R',
        'R': 'Y',
        'S': 'S',
        'W': 'W',
        'K': 'M',
        'M': 'K',
        'B': 'V',
        'D': 'H',
        'H': 'D',
        'V': 'B',
        'N': 'N',
        'a': 't',
        't': 'a',
        'u': 'a',
        'g': 'c',
        'c': 'g',
        'y': 'r',
        'r': 'y',
        's': 's',
        'w': 'w',
        'k': 'm',
        'm': 'k',
        'b': 'v',
        'd': 'h',
        'h': 'd',
        'v': 'b',
        'n': 'n'}
    rcomp= []
    for c in x:
        rcomp.append(rcDict[c])
    rcomp= rcomp[::-1]
    return(''.join(rcomp))

## Open input/output files
## -----------------------
if args.fasta.endswith('gz'):
    fasta= gzip.open(args.fasta)
else:
    fasta= open(args.fasta)

if args.bed.endswith('gz'):
    bed= gzip.open(args.bed)
elif args.bed == '-':
    bed= sys.stdin
else:
    bed= open(args.bed)

if args.fo is None or args.fo == '-':
    fout= sys.stdout
else:
    fout= open(args.fo, 'w')

## Read-in fasta
## -------------
genome= {}
for line in fasta:
    line= line.strip()
    if line.startswith('>'):
        chrom= line.lstrip('>').split()[0]
        if args.verbose:
            sys.stderr.write('Reading: ' + chrom + '\n')
        genome[chrom]= []
        continue
    genome[chrom].append(line)
fasta.close()

for chrom in genome:
    genome[chrom]= ''.join(genome[chrom])

## Extract sequences
## -----------------
for line in bed:
    line= line.strip().split('\t')
    chrom= line[0]
    s= int(line[1])
    e= int(line[2])
    if s >= e:
        sys.stderr.write("Error: Start >= End. Line was: %s\n" %('\t'.join(line)))
        sys.exit(1)
    if args.strand:
        strand= line[5]
    else:
        strand= None

    if args.name:
        try:
            name= line[3]
        except IndexError:
            sys.exit('\nCould not find name field (4th column) in file "%s" for line:\n%s' %(args.bed, '\t'.join(line)))
    else:
        name= chrom + ':' + line[1] + '-' + line[2]
        if args.strand:
            name= name + "(" + strand + ")"

    ## Get sequence    
    faseq= genome[chrom][s:e]
    if args.strand and strand == '-':
        faseq= DNAreverseComplement(faseq)
    
    if args.tab:
        fout.write(name + '\t' + faseq + '\n')
    else:
        fout.write('>' + name + '\n')
        fout.write(faseq + '\n')
    
bed.close()
fout.close()

sys.exit()
