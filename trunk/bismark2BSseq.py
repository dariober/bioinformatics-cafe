#!/usr/bin/env python

import argparse
import sys
import subprocess
import os
import re
import tempfile
import shutil

parser = argparse.ArgumentParser(description="""Convert the methylation report from the Bismark pipeline
to a BSseq obeject.

The input is one or (usually) more files produced by
> bismark_methylation_extractor ... --cytosine_report

Which produces the "genome-wide cytosine methylation output":
# -------------------------------------------------------------------------
SY005_2x2_hmc_Q38_indexed       1       +       0       0       CHH     CTC
SY005_2x2_hmc_Q38_indexed       3       +       1190    16      CHH     CAC
SY005_2x2_hmc_Q38_indexed       5       +       10      1254    CHH     CCC
# -------------------------------------------------------------------------

In general the input files must have in columns
#1: Chrom name
#2: Position
#3: <ignored> Strand
#4: Count metyhlated
#5: Count non-methylated
#6+: Columns ignored

All input files must all be sorted by chrom the by pos and have the same number of rows.
(An error will be issued if inconsistencies are found)

bismark2BSseq.py will extract chrom, pos, #methyalted, #unmethylated from each
input file and produce temporary input files to be read in R and to create a BSseq
object. The output is a R data object containing the BSseq object ready to be read in R.

""")

parser.add_argument('--infile', '-i',           ## Argument name, long and short form
                    required= True,           ## Is this argument required?
                    nargs= '+',
                    help= '''One or more inout files to process. Use - to read the list from stdin
                    ''')  ## Help for this args

parser.add_argument('--outfile', '-o',
                    required= True,
                    help= '''Basename of output file and R object. E.g. -o /my/dir/myBSseq
will create the R data file /my/dir/myBSseq.Rdata which will contain the BSseq object `myBSseq`
                    ''')  ## Help for this args

parser.add_argument('--strip', '-s',
                    required= False,
                    default= None,
                    help= '''A regex passed to re.sub(). This regex will be removed by each
input filename to create sample names. NB: path directories are always stripped.
Default is None which means the entire filename (w/o path) is used. Ignored if --sampleNames
is not None. NB: Sample names must be unqie after having stripped the regex
                    ''')  ## Help for this args

parser.add_argument('--sampleNames', '-n',
                    required= False,
                    default= None,
                    nargs= '+',
                    help= '''List of unique sample names, one for each input file.
Default None (hence use filenames).                    
                    ''')

parser.add_argument('--keeptmp',
                    default= False,
                    action= 'store_true',
                    help= '''Do not delete the temp directory created by python where
the output files are written before reading them into R.
                    ''')


args = parser.parse_args()
# -----------------------------------------------------------------------------

def checkPos(ll):
    """Check chrom and pos in each sublist is the same. I.e. make sure input files
are sorted in the same way and have the same number of rows.
ll is a list of lists each inner list a line of input. checkPos will extract
the first (chrom) and 2nd (pos) item from each inner list and check they are all the same."""
    chromPos= [tuple(x[0:2]) for x in ll]
    if len(set(chromPos)) != 1:
        sys.exit('Files are not ordered in the same way')
        
# -----------------------------------------------------------------------------

## Get input files
if args.infile == ['-']:
    infiles= sys.stdin.readlines()
    infiles= [x.strip() for x in infiles]
else:
    infiles= args.infile
if len(infiles) != len(set(infiles)):
    sys.exit('Duplicate file names in input: %s' %(infiles))

## Create or get sample names
if args.sampleNames is not None:
    sampleNames= args.sampleNames
elif args.strip is not None:
    sampleNames= [re.sub(args.strip, '', os.path.split(x)[1]) for x in infiles]    
else:
    sampleNames= [os.path.split(x)[1] for x in infiles]

if (len(set(sampleNames)) != len(sampleNames)) or (len(set(sampleNames)) != len(infiles)):
    sys.exit('Invalid sample names (not unique and/or not of the same length as input)')

## Prepare output. Each file will be input for BSseq()
tmpdir= tempfile.mkdtemp(prefix= 'bismark2BSseq_')
print('\nTemporary dir set to: %s' %(tmpdir))

covfile= open(os.path.join(tmpdir, 'coverageMatrix.txt'), 'w')
mfile= open(os.path.join(tmpdir, 'MethylationMatrix.txt'), 'w')
chrPosfile= open(os.path.join(tmpdir, 'chromPosTable.txt'), 'w')

namesfile= open(os.path.join(tmpdir, 'nameList.txt'), 'w')
for x in sampleNames:
    namesfile.write(x + '\n')
namesfile.close()

finList= [open(x) for x in infiles]
while True:
    lines= []
    for fin in finList:
        lines.append(fin.readline().strip().split('\t'))
    checkPos(lines)
    if lines[0] == ['']:
        break
    covlist= '\t'.join([str(int(x[3]) + int(x[4])) for x in lines])
    mlist= '\t'.join([str(int(x[3])) for x in lines])
    covfile.write(covlist + '\n')
    mfile.write(mlist + '\n')
    chrPosfile.write('\t'.join([lines[0][0], lines[0][1]]) + '\n')
[x.close() for x in finList]

covfile.close()
mfile.close()
chrPosfile.close()

# ------------------------------------------------------------------------------
# Create R object

rscript= os.path.join(tmpdir, 'Rscript.R')
fin= open(rscript, 'w')
rcmd= """#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(bsseq))
tmpdir<- '%s'
outfile<- '%s'
nameList<- read.table(file.path(tmpdir, 'nameList.txt'), sep= '\t', colClasses= 'character', stringsAsFactors= FALSE)
coverageMatrix<- read.table(file.path(tmpdir, 'coverageMatrix.txt'), sep= '\t', colClasses= 'integer', stringsAsFactors= FALSE)
MethylationMatrix<- read.table(file.path(tmpdir, 'MethylationMatrix.txt'), sep= '\t', colClasses= 'integer', stringsAsFactors= FALSE)
chromPosTable<- read.table(file.path(tmpdir, 'chromPosTable.txt'), sep= '\t', colClasses= c('character', 'integer'), stringsAsFactors= FALSE, col.names= c('chrom', 'pos'))

bsseq<- BSseq(M= as.matrix(MethylationMatrix), Cov= as.matrix(coverageMatrix), pos= chromPosTable$pos, chr= chromPosTable$chrom, sampleNames= as.matrix(nameList)[,1])
save(bsseq, file= outfile)
quit(save= 'no')
""" %(tmpdir, args.outfile)
fin.write(rcmd)
fin.close()
print('\nFinished writing input files. Creating BSseq object...')

cmd= 'Rscript %s' %(rscript)
p= subprocess.Popen(cmd, shell= True, stdout= subprocess.PIPE)
p.wait()
rout = p.stdout.read()
print(rout)

# -------------------------------------------------------------------------------

if not args.keeptmp:
    shutil.rmtree(tmpdir)

sys.exit()

## ls /lustre/sblab/berald01/repository/bismark_out/bismark-*/*.genome-wide_CX_report.txt | bismark2BSseq.py -i - -o bsseq.Rdata -s '\..*'