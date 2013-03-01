#!/usr/local/bin/env python

import sys

docstring='''
DESCRIPTION
    Merge outputs from bismark cytosine report (genome_wide_cytosine_report.pl)
    by summing columns <count methylated>  <count non-methylated> across
    input files.
    Merged file goes to stdout
   
USAGE:
    mergeCytosineReports.py <file1> <file2> ... > merge.report
    ## Read list of files from stdin
    ls *.genome_wide_CpG_report.txt | mergeCytosineReports.py - > merged.report
    
MEMO:
    Input files are expected to have format:
    <chromosome>  <position>  <strand>  <count methylated>  <count non-methylated>  <C-context>  <trinucleotide context>
'''

def mergeLine(lineList):
    """Given a list of lists produce a merged line by summing column 4th and 5th
    Some consistency check also performed.
    Returns: List of strings
    """
    ## Check positions are the same in all files
    pos= [tuple(x[0:3]) for x in lineList]
    if len(set(pos)) != 1:
        sys.exit('Input files do not seem to have the same positions\n%s' %(pos))
    mergeLine= list(pos[0]) 
    methylationSums= [0, 0]
    for line in lineList:
        methylationSums[0] += int(line[3])
        methylationSums[1] += int(line[4])
    mergeLine= mergeLine + methylationSums + lineList[0][5:]
    return([str(x) for x in mergeLine])
# -----------------------------------------------------------------------------

if len(sys.argv) == 1 or sys.argv[1] in ['-h', '--help']:
    sys.exit(docstring)

if sys.argv[1] == '-':
    filelist= sys.stdin.readlines()
    filelist= [x.strip() for x in filelist]
else:
    filelist= sys.argv[1:]

if len(filelist) != len(set(tuple(filelist))):
    sys.exit('Some input files are duplicated.')

fopen= [open(f) for f in filelist]
while True:
    lineList= []
    for fin in fopen:
        line= fin.readline()
        line= line.strip().split('\t')
        lineList.append(line)
    endline= set(tuple([tuple(x) for x in lineList]))
    if endline == set([('',)]):
        break
    else:
        mline= mergeLine(lineList)
        print('\t'.join(mline))
for f in fopen:
    f.close()
sys.exit()

