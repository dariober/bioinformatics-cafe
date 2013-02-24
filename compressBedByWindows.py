#!/usr/bin/env python

import argparse
import subprocess
import tempfile
import sys
import shutil
import os

parser = argparse.ArgumentParser(description= """
DESCRIPTION
    
    Compress a bed file by applying a sliding window with a grouping function
    to the score column (5th col).
    Ouput is sent to stdout and has columns:
    
        <chrom>  <window_start>  <window_end>  <groupedby_score>

EXAMPLE    
    You have a bed file with position of CpGs (one per row), the score column is
    the percentage of Cs methylated. You want to summarize % methylation in windows
    of size 10kb sliding by 1kb:

    compressBedByWindows.py -i methylation.bed -w 10000 -s 1000
    
TODO:
    Use pybedtools, it would allow to read from stdin which is not possible
    at the moment.

""", formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('--input', '-i',
                   required= True,
                   help='''BED file to compress. The 5th is the score column,
i.e. the numeric column to be summarized. input must be SORTED by position, without header.

''')

parser.add_argument('--window_size', '-w',
                   required= True,
                   type= int,
                   help='''Window size to group features. This option passed to
bedtools makewindows.

''')

parser.add_argument('--step_size', '-s',
                   required= False,
                   default= None,
                   type= int,
                   help='''Step size to slide windows. This option passed to
bedtools makewindows. Default step_size= window_size (non-sliding windows)

''')

parser.add_argument('--ops', '-o',
                   required= False,
                   default= 'sum',
                   type= str,
                   help='''Operation to apply to the score column. Default: sum
This option passed to bedtools groupBy.

''')

parser.add_argument('--tmpdir',
                    required= False,
                    default= None,
                   help='''For debugging: Directory where to put the tmp output files.
By default the temp dir is fetched by tempfile.mkdtemp and will be deleted at the end.
With this option the tmp dir will not be deleted. tmpdir will be created if it doesn't
exist.

''')

args = parser.parse_args()

if args.step_size is None:
    args.step_size= args.window_size

if args.tmpdir is None:
    tmpdir= tempfile.mkdtemp(prefix= 'compressBedByWindows_')
else:
    tmpdir= args.tmpdir
    if not os.path.isdir(tmpdir) and not os.path.exists(tmpdir):
        os.makedirs(tmpdir)
    elif not os.path.isdir(tmpdir) and os.path.exists(tmpdir):
        sys.exit('Requested tmp dir %s is a file' %(tmpdir))
    else:
        pass        

basename= os.path.split(args.input)[1]
        
cmd= """
## 1. Get extremes of each chrom
groupBy -g 1 -c 2,3 -o min,max -i %(inputBed)s > %(tmpdir)s/%(basename)s.chromSize &&

## 2. Divide each chrom in wndows
bedtools makewindows -b %(tmpdir)s/%(basename)s.chromSize -w %(w)s -s %(s)s > %(tmpdir)s/%(basename)s.windows &&

## 3. Assign bed features to windows
intersectBed -a %(tmpdir)s/%(basename)s.windows -b %(inputBed)s -wa -wb > %(tmpdir)s/%(basename)s.intsct && 

## 4. Summarize windows
groupBy -g 1,2,3 -c 8 -o %(ops)s -i %(tmpdir)s/%(basename)s.intsct
""" %{'tmpdir': tmpdir, 'inputBed': args.input, 'basename': basename, 'w':args.window_size, 's': args.step_size, 'ops': args.ops}

fsh= open(os.path.join(tmpdir, basename) + '.sh', 'w')
fsh.write(cmd)
fsh.close()

p= subprocess.Popen(cmd, shell= True)
p.communicate()

if args.tmpdir is None:
    shutil.rmtree(tmpdir)
sys.exit()