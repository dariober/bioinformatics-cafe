#!/usr/bin/env python

"""
Note: f.tell() and f.seek() on gzip files is extremely slow. Do not use
gzip files.
"""

import argparse
import os
import sys
import glob
import re
import subprocess
from itertools import imap
from operator import itemgetter
import heapq

parser = argparse.ArgumentParser(description= """
DESCRIPTION
    Produce union loci and count matrices from one or more files from
    mpileup2methylation.py.

EXAMPLE USAGE
    BSmatrix.py -i *.mcall.bdg -p bsmat -o ./bsdata -s '\.*'

NOTES:
    Glob characters in input are expanded and duplicate files are included only
    once in output. The order of the files in output is alphanumeric based on the
    file name (only name w/o directory path).

CAVEATS & ASSUMPTIONS:
    - Input files are sorted by chromosome and position. There is a check for
    correct chrom order but not for position
    For sorting by chromosome and position use:
    sort -k1,1 -k2,2n -k3,3n
    For sorting chromosome only (positions already sorted):
    sort -k1,1 -s
    - Input must be unzipped (processing gzip possible but slooooower)

""", prog= 'BSmatrix.py', formatter_class= argparse.RawDescriptionHelpFormatter) # , formatter_class= argparse.RawTextHelpFormatter

# -----------------------------------------------------------------------------

parser.add_argument('--version', action='version', version='%(prog)s 0.1.0a')

parser.add_argument('--input', '-i',
                   required= True,
                   default= [],
                   nargs= '+',
                   help='''List of input files to cross-tabulate to matrices.
Glob (wildcard) patterns are expanded.
                   ''')

parser.add_argument('--prefix', '-p',
                   required= True,
                   help='''String of prefix to use for the output files.
                   ''')

parser.add_argument('--outdir', '-d',
                   required= False,
                   default= '.',
                   help='''Output directory default to current dir. Created if it doesn't
exist.
                   ''')

parser.add_argument('--strip', '-s',
                   required= False,
                   default= None,
                   help='''A regular expression to strip from the input file names.
The resulting strings will be used as names in output.
                   ''')
# -----------------------------------------------------------------------------

def listOpener(inputList):
    """Return a list of open file handles for reading from the list of input file
    names.
    """
    input_fin= []
    for x in inputList:
        try:
            fin= open(x)
#            if x.endswith('.gz'):
#                fin= gzip.open(x)
#            else:
#                fin= open(x)
        except IOError:
            sys.exit('Cannot open file "%s" for reading.' %(x))
        input_fin.append(fin)
    return(input_fin)
    
def extract_key(line):
    """Extract keys from each input line (bdg file)
    """
    line= line.split('\t')
    key= [line[0], int(line[1]), int(line[2])]
    return(key)

def bdgLine_to_locusLine(line):
    """Convert a line from a bdg file into a line formatted for locus map.
    Input and output are is a *strings*
    """
    line= line.strip().split('\t')
    xline= line[0:3] + ['_'.join(line[0:3]), '.', '+']
    xline= '\t'.join(xline)
    return(xline)
# -----------------------------------------------------------------------------
args = parser.parse_args()

## Check and open INPUT
## --------------------
input= []
for x in args.input:
    xg= glob.glob(x)
    input += xg
input= sorted(list(set(input)))
if input == []:
    sys.exit('No files in input')

input_fin= listOpener(input)

## Prepare LIBRARY NAMES
## ---------------------

outnames= []
for x in input:
    fname= os.path.split(x)[1]
    if args.strip:
        outnames.append(re.sub(args.strip, '', fname))
    else:
        outnames.append(fname)
if len(outnames) != len(set(outnames)):
    sys.exit("Duplicate output IDs found in %s" %(outfiles))

print('\nInput file => library name:')
for i,o in zip(input, outnames):
    print(i + ' =>  ' + o)

## Prepare OUTPUT FILES
## --------------------
if not os.path.exists(args.outdir):
    try:
        os.makedirs(args.outdir)
    except OSError:
        sys.exit('Unable to create output directory "%s"' %(args.outdir))

if '/' in args.prefix or args.prefix.strip().startswith('.') or args.prefix.startswith(' ') or args.prefix.strip() == '':
    sys.exit('Invalid prefix. Got "%s" ' %(args.prefix))
    
outprefix= os.path.join(args.outdir, args.prefix)

floci= outprefix + '.loci'
fpctmet= outprefix + '.pct_met.mat'
fcntmet= outprefix + '.cnt_met.mat'
ftotreads= outprefix + '.tot_reads.mat'

try:
    fout_loci= open(floci, 'w')
    fout_pct_met= open(fpctmet, 'w')
    fout_cnt_met= open(fcntmet, 'w')
    fout_tot_reads= open(ftotreads, 'w')
except IOError:
    sys.exit('I cannot create output files in dir "%s" with prefix "%s".' %(args.outdir, args.prefix))

fout_loci.write('\t'.join(['chrom', 'start', 'end', 'locus', 'score', 'strand']) + '\n')

fout_pct_met.write('\t'.join(['locus'] + outnames) + '\n')
fout_cnt_met.write('\t'.join(['locus'] + outnames) + '\n')
fout_tot_reads.write('\t'.join(['locus'] + outnames) + '\n')

print('\nOutput files:\n%s' %('\n'.join([floci, fpctmet, fcntmet, ftotreads])))
## ----------------------------
## 1. Create union of positions
## ----------------------------

## See http://stackoverflow.com/questions/12460943/merging-pre-sorted-files-without-reading-everything-into-memory
print('\nGenerating union of all positions...')

#with input_fin[0] as f0, input_fin[1] as f1, input_fin[2] as f2, input_fin[3] as f3, input_fin[4] as f4, input_fin[5] as f5:
#with input_fin as sources:
#with fout_loci as dest:
decorated = [
    ( ( extract_key(line), bdgLine_to_locusLine(line) ) for line in f if line.strip() != '' )
    for f in input_fin]
merged = heapq.merge(*decorated)
undecorated = imap(itemgetter(-1), merged)
thisline= ''
print('Outputting')
for x in undecorated:
    if x != thisline:
        fout_loci.write(x + '\n')
        nlines += 1
        if nlines % 500000 == 0:
            print('%s lines to locus file' %(nlines))
#        if nlines > 50000:
#            break
    thisline= x
fout_loci.close()
for x in input_fin:
    x.close()

## 2. Re-pass through all the files to generate matrices (rows are from the union)
## -----------------------------------------------------------------------------
print('\nGenerating matrices...')
fout_loci= open(floci)
x= fout_loci.readline()
input_fin= listOpener(input)

chromList= [None] * len(input_fin) ## This is to check that current chrom is >= previous
m= 0
for line in fout_loci:
    line= line.strip().split('\t')
    pos= [line[0], int(line[1]), int(line[2])]
    locusName= line[3]
    totLine= [locusName]
    cntmetLine= [locusName]
    pctmetLine= [locusName]
    for i in range(0, len(input_fin)):
        fin= input_fin[i]
        t= fin.tell()
        xline= fin.readline().strip().split('\t')        
        if xline != ['']:
            xpos= [xline[0], int(xline[1]), int(xline[2])]
            if xpos == pos:
                pctmetLine.append(xline[3])
                cntmetLine.append(xline[4])
                totLine.append(xline[5])
            elif xpos > pos:
                """There is no count at this position, put 0 and hold back the file
                reading"""
                fin.seek(t)    
                pctmetLine.append('0')
                cntmetLine.append('0')
                totLine.append('0')
            else:
                sys.exit('Something wrong')
            ## Check chroms are in the right order
            if chromList[i] < xline[0]:
                chromList[i]= xline[0]
            elif chromList[i] == xline[0]:
                pass
            elif chromList[i] > xline[0]:
                sys.exit('\n\nChromosomes are not in the right order in file "%s": %s (current) < %s (previous)\nHint: sort with `sort -k1,1 -k2,2n -k3,3n` or `sort -k1,1 -s`\n' %(input[i], xline[0], chromList[i]))
            else:
                sys.exit('Unexpected case')
        else:
            pctmetLine.append('0')
            cntmetLine.append('0')
            totLine.append('0')
    m += 1
    if m % 1000000 == 0:
        print('%s/%s lines processed' %(m,nlines))
    fout_pct_met.write('\t'.join([str(x) for x in pctmetLine]) + '\n')
    fout_cnt_met.write('\t'.join([str(x) for x in cntmetLine]) + '\n')
    fout_tot_reads.write('\t'.join([str(x) for x in totLine]) + '\n')

for fin in input_fin:
    xline= fin.readline().strip()
    if xline != '':
        print('Not all lines read from one or more input files! New line read:')
        print(xline)
        sys.exit(1)
    fin.close()
fout_loci.close()
fout_cnt_met.close()
fout_pct_met.close()
fout_tot_reads.close()

sys.exit()
