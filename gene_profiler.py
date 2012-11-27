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
    Divides each input bed feature in n equally sized windows and counts the overlaps
    with the input bam. Bins are then grouped to produce a coverage profile.

OUTPUT:
    Output is a tab separated with header line:
    bin: Bin number (see below)
    mean_rpkm: Mean RPKM for this bin number (see below).
    mean_count: Mean of row count (sum_count/no_bins)
    sum_count: Sum of row counts
    mean_size: Mean size of this bin number
    no_bins: Number of bins (i.e. number of bed features *after* merging)

NOTES:
    - Bins are numbered so that bin number 1 is always the most upstream and bin
        n the most dowmnstream. E.g.
        If a gene is on - strand, bin 1 is on TSS and will be the rightmost.
    - Input BED file must be a BED6 (additional columns ignorerd).
    - Rpkm is calculated using the total number of reads mapping to bed features
        (not the total number of reads in the BAM)
    - Before counting, bed features are merged in a strand specific way. Threfore overalapping
        features on opposite strands are double counted!
    
REQUIRES:
    bedtools 
    bed_windows.py (http://code.google.com/p/bioinformatics-misc/source/browse/trunk/bed_windows.py)

EXAMPLES:

TODO:
    Extend rads by x bases using slopBed (see snippet in the code).
    Exit with error if any of the shell commands don't return clean.
    
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
                    help="""Name for the output file. 
                    """)

parser.add_argument('-w', '--nwins',
                    type= int,
                    default= 100,
                    help="""Divide each bed feature in this many windows. Default 100. 
                    """)
parser.add_argument('--skip',
                    type= int,
                    default= 0,
                    help="""Skip these many lines from input bed. E.g. set 1 for
skipping the header. Default 0.
                    """)

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

## Awk will make bed6 format
print('\n\nMerging overalpping features:\n')
cmd= '''tail -n+%s %s | sort -k1,1 -k2,2n | mergeBed -s -i stdin | awk '{print $1"\\t"$2"\\t"$3"\\tmerged\\t"$3-$2"\\t"$4}' > %s''' %(args.skip, args.bed, mergebed) 
print(cmd)
p= subprocess.Popen(cmd, shell= True); p.wait()

# ------------------[ Prepare windows ]-----------------------------------------

print('\n\nPreparing windows\n')
## To comply with BED6, awk will place the bin number on column 4 and the strand on column 6. Column 5 will be window size
cmd= '''bed_windows.py --nwinds %(nwins)s -r %(bedfeatures)s | awk '{print $1"\\t"$2"\\t"$3"\\t"$9"\\t"$3-$2"\\t"$8}' > %(winds)s ''' %{'bedfeatures':mergebed, 'winds': winds, 'nwins': args.nwins}
print(cmd)
p= subprocess.Popen(cmd, shell= True); p.wait()

""" This part uses bedtools to make windows
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
"""
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

## MEMO: File geneprofilerpkm has format:
## 1) Chrom
## 2) Bin start
## 3) Bin end
## 4) Bin number
## 5) Bin size
## 6) Strand
## 7) Bin RPKM
## 8) The number of features in A that overlapped the B interval.
## 9) The number of bases in B that had non-zero coverage.
## 10) The length of the entry in B.
## 11) The fraction of bases in B that had non-zero coverage.

# -------------------------[ GroupBy (summarize) ]-----------------------------
print("\n\nSummarize coverage\n")
header= '\t'.join(['bin', 'mean_rpkm', 'mean_count', 'sum_count', 'mean_size', 'no_bins'])
fout= open(args.out, 'w')
fout.write(header + '\n')
fout.close()
cmd= '''bedtools groupby -i %(geneprofile)s -g 4 -c 7,8,8,10,10 -o mean,mean,sum,mean,count  >> %(profileGroupby)s''' %{'geneprofile': geneprofilerpkm, 'profileGroupby':args.out}
print(cmd)
p= subprocess.Popen(cmd, shell= True); p.wait()

# -------------------------------[ Clean-up ]----------------------------------
if not args.keeptmp:
    shutil.rmtree(tmpdir)
else:
    print('\n\nTemp files have not been removed and are in:\n%s/\n' %(tmpdir))
sys.exit()