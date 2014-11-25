#!/usr/bin/env python

import argparse
import subprocess
import tempfile
import sys
import os
import shutil
import scipy.stats
import numpy
import atexit

parser = argparse.ArgumentParser(description= """
DESCRIPTION
    Compute the read enrichment in target intervals relative to local background.
    
    Typical use case: A ChIP-Seq experiment on a sample returns a number of regions
    of enrichemnt. We want to know how enriched these regions are in a different
    sample. Note that enrichment is quantified relative to the local background
    not relative to an input control.
    
OUTPUT:
    bed file with header and columns:
1. chrom
2. start
3. end
4. targetID
5. flank_cnt
6. target_cnt
7. flank_len
8. target_len
9. log10_pval
10. log2fc

EXAMPLE
    localEnrichmentBed.py -b rhh047.bam -t rhh047.macs_peaks.bed -g genome.fa.fai -bl blacklist.bed > out.bed
    
    ## Using pipes:
    samtools view -u rhh047.bam chr18 \
    | localEnrichmentBed.py -b - -t rhh047.macs_peaks.bed -g genome.fa.fai > out.bed

    
Useful tip: Get genome file from bam file:

    samtools view -H rhh047.bam \
    | grep -P "@SQ\tSN:" \
    | sed 's/@SQ\tSN://' \
    | sed 's/\tLN:/\t/' > genome.txt
    
REQUIRES:
    - bedtools suite
    - numpy, scipy 
    
""", formatter_class= argparse.RawTextHelpFormatter, prog= os.path.basename(__file__))

parser.add_argument('--target', '-t',
                   required= True,
                   help='''Target bed file where enrichment is to be computed.
Use - to read from stdin.
                ''')

parser.add_argument('--bam', '-b',
                   required= True,
                   help='''Bam file of the library for which enrichment is to be
computed. Use - to read from stdin.
                   ''')

parser.add_argument('--genome', '-g',
                   required= True,
                   help='''A genome file giving the length of the chromosomes.
A tab separated file with columns <chrom> <chrom lenght>.
NB: It can be created from the header of the bam file (see tip above).
                   ''')

parser.add_argument('--slop', '-S',
                   default= '5.0',
                   help='''Opt passed to slopBed to define the surrounding of the
target regions:
If `int` each target will be extended left and right this many bases.
If `float` each target is extended left and right x times its size.
E.g. 5.0 (default) extends each target regions 5 times its length.
                   ''')

parser.add_argument('--blacklist', '-bl',
                   required= False,
                   help='''An optianl bed file of regions to ignore to compute
the local background. These might be unmappable regions with 0-counts which would
inflate the target enrichment.
                   ''')

parser.add_argument('--tmpdir',
                   required= False,
                   help='''Temp dir to use for intermediate files. If not set
python will get one. A subdir will be created here.
                   ''')

parser.add_argument('--keeptmp',
                   action= 'store_true',
                   help='''If set, the tmp dir is not deleted at the end of the
job (useful for debugging).
                   ''')

parser.add_argument('--verbose', '-V',
                   action= 'store_true',
                   help='''Print to stderr the commands that are executed. 
                   ''')

parser.add_argument('--version', action='version', version='%(prog)s 0.1')

args= parser.parse_args()
# ------------------------------------------------------------------------------

def prepareTargetBed(inbed, outbed, verbose= False):
    """Create a bed file where the fourth column is a unique identifier and the
    5th column is marked 'target'.
    inbed:
        Name of input bed file
    outbed:
        Name for output bed file        
    
    Returns True on exit.
    """
    if verbose:
        sys.stderr.write('\nWriting target bed to "%s"\n' %(outbed))
        
    if inbed == '-':
        fin= sys.stdin
    else:
        fin= open(inbed)
    id= 1
    outtarget= open(outbed, 'w')
    for line in fin:
        line= line.strip().split('\t')
        nr= str(id) + '\ttarget'
        outtarget.write('\t'.join(line[0:3] + [nr]) + '\n')
        id += 1
    outtarget.close()
    fin.close()
    return(True)

def prepareFlankingRegions(targetBed, slop, genome, blacklistBed= None, flankingBed= 'flankingBed', verbose= False):
    """Extend the target regions to define the flanking regions to use as
    background. If a flanking region overlaps other targets or a blacklisted
    region, the intersection is removed (so flanking regions are not necessarily
    one left and one right).
        
    targetBed:
        Name of file of target regions to extend.
        Typically from prepareTargetBed
    slop:
        int or float to extend the targets. If int extend by this many bp. If
        float, extend the target by this many times its length.
    genome:
        Name of genome file passed to slopBed
    flankingBed:
        Name of output file for flanking regions.
    blacklistBed:
        Name of file of blavklist regions.
    """
    ## Determine if "slop" is int or float:
    try:
        S= int(slop)
    except ValueError:
        try:
            S= float(slop)
        except:
            sys.exit('\nInvalid slop. Must be int or float. Got: "%s".\n' %(slop))
    if type(S) == int:
        pct= ''
    else:
        pct= '-pct'
    
    # Prepare subcommand for blacklist. If BL is not given, do nothing
    if blacklistBed:
        cmdBL= """| subtractBed -a - -b %(blacklist)s \\""" %{'blacklist': blacklistBed}
    else:
        cmdBL= "\\"
    
    # Flanking
    cmd= """
slopBed %(pct)s -l %(S)s -r %(S)s -g %(genome)s -i %(target)s \\
| subtractBed -a - -b %(target)s \\
%(cmdBL)s
| awk 'BEGIN{OFS="\\t"}{print $1, $2, $3, $4, "flank"}' > %(flankingBed)s
""" %{'pct':pct, 'S':S, 'genome': genome, 'target': targetBed, 'cmdBL': cmdBL, 'flankingBed': flankingBed}
    
    if verbose:
        sys.stderr.write(cmd)
    
    p= subprocess.Popen(cmd, shell= True, stdout= subprocess.PIPE, stderr= subprocess.PIPE)
    stdout, stderr= p.communicate()
    if p.returncode != 0:
        print(stderr)
        print('Exit code %s' %(p.returncode))
        sys.exit(1)
    return(True)

def countReadsInInterval(inbam, beds, countTable, tmpdir, verbose= False):
    """Count reads overalapping intervals in bed file and group by target id
    (4th columns) and region type (5th column).
    inbam:
        Name of bam file to get reads from
    beds:
        List of bed file names where reads will be counted. They will be concatenated.
        Typically these are targetBed and flankingBed.
    countTable:
        Name of output file where counts will be stored.
    """
    
    cmd= """
cat %(beds)s > %(catbeds)s &&
coverageBed -counts -abam %(bam)s -b %(catbeds)s \\
| awk 'BEGIN{OFS="\\t"}{print $0, $3-$2}' \\
| sort -s -k4,4n -k5,5 \\
| groupBy -g 4,5 -c 6,7 -o sum,sum -prec 12 > %(countTable)s
""" %{'beds': ' '.join(beds), 'catbeds': os.path.join(tmpdir, 'cat.bed'), 'bam': inbam, 'countTable':countTable}
    
    if verbose:
        sys.stderr.write(cmd)

    p= subprocess.Popen(cmd, shell= True, stdout= subprocess.PIPE, stderr= subprocess.PIPE, stdin= subprocess.PIPE)
    if inbam == '-':
        for line in sys.stdin:
            p.stdin.write(line)
    
    stdout, stderr= p.communicate()
    if p.returncode != 0:
        print(stderr)
        print('Exit code %s' %(p.returncode))
        sys.exit(1)

    return(True)

def countsToDict(fin):
    """Read a set of line from file connection fin and return a tuple as
    (int, {'flank': {'cnt': int, 'len': int}, 'target': {'cnt': int, 'len': int}})
    """
    locus= {}
    while True:
        x= fin.tell()
        line= fin.readline()
        if line != '':
            line= line.strip().split('\t')
            thisLocus= line[0]
            if locus == {} or thisLocus == current:
                current= line[0]
                locus[line[1]]= {'cnt': int(line[2]), 'len': int(line[3])}
            else:
                fin.seek(x)
                return((int(current), locus))
        elif locus != {}:
            return((int(current), locus))
        else:
            return(None)

def localEnrichment(countTuple):
    """Take a tuple composed of (ID, Dict) to extract chi^2 stats for enrichment
    and log2 fold change.
    countTuple:
        Tuple with ID and dictionary to extract counts. Typically returned by
        countsToDict.
    Return: Tuple with dictionary updated with -log10(pvalue) and log2 FC
    """
    # ct= {'flank': {'cnt': 142, 'len': 7110}, 'target': {'cnt': 111, 'len': 711}}
    ct= countTuple[1]
    cnt= [ct['flank']['cnt'], ct['target']['cnt']]
    length= [ct['flank']['len'], ct['target']['len']]
    
    try:
        ct['log10_pval']= -numpy.log10(scipy.stats.chi2_contingency([cnt, length])[1])
    except ValueError:
        ct['log10_pval']= None
        
    if ct['flank']['cnt'] == 0 or ct['flank']['len'] == 0:
        ct['log2fc']= None
    else:
        ct['log2fc']= numpy.log2(
             float((ct['target']['cnt']) / float(ct['target']['len'])) /
            (float(ct['flank']['cnt']) / float(ct['flank']['len'])))
        
    countTupleUp= (countTuple[0], ct)
    return(countTupleUp)

# ------------------------------------------------------------------------------
## Get tmp dir:
tmpdir= tempfile.mkdtemp(prefix= 'localEnrichmentBed_', dir= args.tmpdir)
if not args.keeptmp:
    atexit.register(shutil.rmtree, tmpdir)

if args.verbose:
    sys.stderr.write('\nWorking tmp dir: "%s"\n' %(tmpdir))

targetBed= os.path.join(tmpdir, 'target.bed')
prepareTargetBed(args.target, targetBed, verbose= args.verbose)

flankingBed= os.path.join(tmpdir, 'flanking.bed')
prepareFlankingRegions(targetBed= targetBed, slop= args.slop, genome= args.genome,
    blacklistBed= args.blacklist, flankingBed= flankingBed, verbose= args.verbose)

countTable= os.path.join(tmpdir, 'countTable.txt')
countReadsInInterval(inbam= args.bam, beds= [targetBed, flankingBed],
    countTable= countTable, tmpdir= tmpdir, verbose= args.verbose)

header= 'chrom, start, end, targetID, flank_cnt, target_cnt, flank_len, target_len, log10_pval, log2fc'.replace(', ', '\t')
print(header)
fcount= open(countTable)
ftarget= open(targetBed)
while True:
    cntDict= countsToDict(fcount)
    if not cntDict:
        break
    id, ct= localEnrichment(cntDict)
    outline= ftarget.readline().strip().split('\t')[0:4]
    if str(outline[3]) != str(id):
        sys.exit('Target file and count table not in sync! Got IDs: "%s" and "%s"' %(outline[3], id))
    outline.append(str(ct['flank']['cnt']))
    outline.append(str(ct['target']['cnt']))
    outline.append(str(ct['flank']['len']))
    outline.append(str(ct['target']['len']))
    outline.append(str(ct['log10_pval']))
    outline.append(str(ct['log2fc']))
    print('\t'.join(outline))

fcount.close()
ftarget.close()

sys.exit()
