#!/usr/bin/env python

import os
import argparse
import sys
import re
import tempfile
import shutil
import subprocess

parser = argparse.ArgumentParser(description= """

DESCRIPTION:

    
EXAMPLES:
    

TODO:
      
""", formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('--bed',
                    type= str,
                    required= True,
                    help="""Input bed files for which the profile should be
computed.
                    """)

parser.add_argument('--bam',
                    type= str,
                    required= True,
                    help="""Input bam file.
                    """)

parser.add_argument('-o', '--out',
                    type= str,
                    required= True,
                    help="""Name of the output file. 
                    """)

parser.add_argument('-w', '--nwins',
                    type= int,
                    default= 100,
                    help="""Divide each bed feature in this many windows. Default 100. 
                    """)

#parser.add_argument('--groupbyOps',
#                    type= str,
#                    default= 'mean,sum',
#                    help="""The operation that will be applied to the each window.
#This arg passed to bedtools "groupBy -ops". See "groupBy -h" for possible options.
#Default: mean
#                    """)

parser.add_argument('--tmpdir',
                    type= str,
                    default= None,
                    required= False,
                    help="""Temp directory where to put intermediate files. It
will be created if it doesn't exist. Default: Python will found one using tempfile.mkdtemp
                    """)

parser.add_argument('--keeptmp',
                    action= 'store_true',
                    default= False,
                    help="""Should the temporary directory be removed? With this flag
it will not be removed.
                    """)

args = parser.parse_args()

# -----------------------------------------------------------------------------

def totReadsCoverageBed(covbed):
    """Get the total number of reads mapped to the windows. This is just the
    sum of column 7 from the output of coverageBed.
    ARGS:
    covbed:
        The output of coverageBed. The column to sum is expected to be the 7th
    RETURN:
        int of sum of counts
    """
    fin= open(covbed)
    totCount= 0
    for line in fin:
        totCount += int(line.strip().split('\t')[6])
    fin.close()
    return(totCount)

def rpkmCoveragebed(covbed, N):
    """Normalize by RPKM the count in each bed feature produced by coverageBed.
    The column to normalize is expected to be the 7th.
    ARGS:
    covbed:
        The output of coverageBed. The column to sum is expected to be the 7th
    N:
        Normalize by this total count. Typically this is the total number of features
        mapped to the features (output of totReadsCoverageBed).
    
    RETURN:
        A file named <covbed>.rpkm identical to covbed but with the normalized count at position 7 and the name of the file itself     
        """
    fin= open(covbed)
    outfile= covbed + '.rpkm'
    fout= open(outfile, 'w')
    for line in fin:
        line= line.strip().split('\t')
        L= int(line[2]) - int(line[1])
        C= float(line[6])
        rpkm= (1000000000.0 * C / (N * L))
        line.insert(6, str(rpkm))
        fout.write('\t'.join(line) +'\n')
    fin.close()
    fout.close()
    return(outfile)
    
# -----------------------------------------------------------------------------

# Prepare output
if args.tmpdir is None:
    tmpdir= tempfile.mkdtemp(prefix= 'tmp_gene_profiler_')
else:
    tmpdir= args.tmpdir
    if not os.path.exists(tmpdir):
        os.makedirs(tmpdir)
    
mergebed= os.path.join(tmpdir, args.bed + '.merged.bed')
windplus= os.path.join(tmpdir, 'windows.plus.bed') ## $(mktemp tmp/${geneprofile}.windows.plus.XXXX.bed)
windminus= os.path.join(tmpdir, 'windows.minus.bed')
winds= os.path.join(tmpdir, 'windows.bed')  ## $(mktemp tmp/${geneprofile}.windows.XXXX.bed)
geneprofile= os.path.join(tmpdir, 'geneprofile.bed')

# ----------------------[ Merge overalapping features ]------------------------

cmd= 'sort -k1,1 -k2,2n %s | mergeBed -s -i stdin > %s' %(args.bed, mergebed) 
print(cmd)
p= subprocess.Popen(cmd, shell= True); p.wait()

# ------------------[ Prepare windows ]-----------------------------------------

print("\n\nProduce windows for plus and minus strand\n")
cmd= '''grep -P '\\t\+$' %(bedfeatures)s | cut -f 1-3 | sort -k1,1 -k2,2n -k3,3n | \
     uniq | bedtools makewindows -b stdin -n %(nwins)s -i winnum | awk '{print $1"\\t"$2"\\t"$3"\\t"$4"\\t.\\t+"}' > %(windplus)s ''' %{'bedfeatures':mergebed, 'windplus': windplus, 'nwins': args.nwins} ## Bin 1 is on the TSS
print(cmd)
p= subprocess.Popen(cmd, shell= True); p.wait()
cmd= '''grep -P '\\t\-$' %(bedfeatures)s | cut -f 1-3 | sort -k1,1 -k2,2n -k3,3n | \
     uniq | bedtools makewindows -b stdin -n %(nwins)s -i winnum | awk -v N=$%(nwins)s '{print $1"\\t"$2"\\t"$3"\\t"(N-$4)+1"\\t.\\t-"}' > %(windminus)s ''' %{'bedfeatures':mergebed, 'windminus': windminus, 'nwins': args.nwins}
print(cmd)
p= subprocess.Popen(cmd, shell= True); p.wait()
cmd= '''cat %(windplus)s %(windminus)s > %(winds)s''' %{'windplus':windplus, 'windminus':windminus, 'winds': winds}
print(cmd)
p= subprocess.Popen(cmd, shell= True); p.wait()

# ------------------[ Produce coverage ]----------------------------------------

print("\n\nProfile coverage\n")

## Extend reads by X many bases. See also http://gabriele-bucci.blogspot.co.uk/2012/02/create-bigwig-file-from-bam.html
## samtools view -b %(bam)s | bamToBed | slopBed -r 200 -s -i stdin -g chrSizes.genome | coverageBed -a stdin -b %(winds)s

cmd= '''coverageBed -abam %(bam)s -b %(winds)s | sort -S 20%% -k4,4n > %(geneprofile)s''' %{'bam':args.bam, 'winds': winds, 'geneprofile': geneprofile}
print(cmd)
p= subprocess.Popen(cmd, shell= True); p.wait()

print("\n\nNormalizing read counts...")
totreads= totReadsCoverageBed(geneprofile)
geneprofilerpkm= rpkmCoveragebed(geneprofile, totreads)
print("Reads mapped to windows: %s\nNormalized read count in %s" %(totreads, geneprofilerpkm))

## MEMO: Default output of coverageBed (2.17) has the input bed with last 4 columns being:
#          1) The number of features in A that overlapped the B interval.
#	   2) The number of bases in B that had non-zero coverage.
#	   3) The length of the entry in B.
#	   4) The fraction of bases in B that had non-zero coverage.


print("\n\nSummarize coverage\n")
cmd= '''bedtools groupby -i %(geneprofile)s -g 4 -c 8,8,10,10 -o mean,sum,mean,count  > %(profileGroupby)s''' %{'geneprofile': geneprofilerpkm, 'profileGroupby':args.out}
print(cmd)
p= subprocess.Popen(cmd, shell= True); p.wait()

if not args.keeptmp:
    shutil.rmtree(tmpdir)
else:
    print('\n\nTemp files have not been removed and are in:\n%s/\n' %(tmpdir))
sys.exit()