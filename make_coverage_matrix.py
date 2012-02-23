#! /usr/local/bin/python

import argparse
import sys
import subprocess
import os
import re
# -------------------------[ Parse argumets ]----------------------------------

parser = argparse.ArgumentParser(description="""Produce a matrix of number of reads per feature (peak). Input is all BED files in indir produced by coverageBed (usually the 4th column from right). Feature (peak) ID is in column 4. NB: Input files should have the same number of rows and sorted in the same way. """)
parser.add_argument('--indir', '-i',           ## Argument name, long and short form
                    default= '.',              ## Default to current working dir
                    required= False,           ## Is this argument required? 
                    help= 'Input directory where the BED files to process are. Default to current dir.')  ## Help for this args
parser.add_argument('--index', '-x',           ## Argument name, long and short form
                    required= False,           ## Is this argument required?
                    type= int,
                    default= -4,
                    help= 'The index position of the column to extract. 0 based. E.g. -x 1 extracts 2nd column, -x -2 extracts 2nd column from right (=2nd last). Default is suitable for coverageBed ran with default args.')  ## Help for this args
parser.add_argument('--extension', '-e',           ## Argument name, long and short form
                    required= False,           ## Is this argument required?
                    type= str,
                    default= '\.bed$',
                    help= 'Parse only the file in indir matching this regular expression')  ## Help for this args
parser.add_argument('--subregex', '-s',           ## Argument name, long and short form
                    required= False,           ## Is this argument required?
                    type= str,
                    default= None,
                    help= 'Pass this regex to re.sub() to remove it from the file name')  ## Help for this args
parser.add_argument('--rowheaders', '-r',
                    required= False,
                    type= int,
                    nargs='+',
                    default= [3],
                    help= 'One or more column indexes to use as row headers. Default is 3 (4th column of bed file that is feature name)')  ## Help for this args

args = parser.parse_args()

# ---------------------------[ Define functions ]------------------------------

def get_pathbed(beddir):
    """
    Get all the bed files in directory beddir and return the files as list with beddir path included
    E.g.
    get_pathbed(solexa_pipeline)
    >>> ['solexa_pipeline/myfile1.bed', 'solexa_pipeline/myfile2.bed', ...]
    """
    allfiles= os.listdir(beddir)
    allfiles.sort()
    bedfiles= []
    for f in allfiles:
        if re.search(args.extension, f):
            f= os.path.join(beddir, f)
            bedfiles.append(f)
    return(bedfiles)


# ------------------------------------------------------------------------------

bedfiles= get_pathbed(args.indir)
if len(bedfiles) == 0:
    sys.exit()
bedfiles.sort()

infiles=[]
for i in range(0, len(bedfiles)):
    infile= 'infile'+ str(i)
    vars()[infile]= open(bedfiles[i])
    infiles.append(vars()[infile])

filenames= [os.path.split(x)[1] for x in bedfiles]
if args.subregex is not None:
    filenames_header= [re.sub(args.subregex, '', x) for x in filenames]
else:
    filenames_header= filenames
preheader= ['col_id_' + str(i) for i in range(1, len(args.rowheaders)+1)]
header_row= preheader + filenames_header
header_row= '\t'.join(header_row)
print(header_row)

done= 0
while True:
    matline= []
    for i in range(0, len(infiles)):
        if i == 0:
            " Extract both peak name and count "
            line= (infiles[i]).readline().strip()
            if line == '':
                done= 1
                break
            line= line.strip().split('\t')
            rowheaders= [line[int(x)] for x in args.rowheaders]
            line= rowheaders + [line[args.index]]
            line= '\t'.join(line)
            matline.append(line)
        else:
            " Extract only the count "
            line= (infiles[i]).readline().strip()
            if line == '':
                done= 1
                break
            line= line.strip().split('\t')
            line= line[args.index]
            matline.append(line)
    matline= '\t'.join(matline)
    if done == 1:
        break
    print(matline)

for f in infiles:
    f.close()

sys.exit()
