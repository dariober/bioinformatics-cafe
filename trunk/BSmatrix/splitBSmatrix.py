#!/usr/bin/env python

docstring= """Split output from BSmatrix to individual chromosmes
IMPORTANT
    chrom is either the first column (.loci) or the first part of the locus
    name after the underscore (e.g.: chr10_1_2 => chr10)

USAGE:
    splitBSmatrix.py <basename> <outdir> 

basename:
    Basename to recognize matrix files (e.g. "genome" will look for "genome.loci",
    "genome.tot_reads.mat", etc)
outdir:
    Output directory defualt pwd. One subdir per chromosme will be created. Files
    under subdirs will have the same prefix as basename.
"""

import sys
import os

if len(sys.argv) == 1 or len(sys.argv) > 3 or sys.argv[1] in ['-h', '--help']:
    print(docstring)
    sys.exit(1)

basename= sys.argv[1]
if len(sys.argv) == 3:
    outdir= sys.argv[2]
else:
    outdir= os.getcwd()

infiles= [basename + '.loci', basename + '.pct_met.mat', basename + '.cnt_met.mat', basename + '.tot_reads.mat']
fin= []
for f in infiles:
    print('Opening %s' %(f))
    try:
        fin.append(open(f))        
    except:
        sys.exit('Cannot open input file %s' %(f))

headers= []
for f in fin:
    header= f.readline().strip()
    headers.append(header)

prevChrom= None
ndirs= 1
for i in range(0, len(fin)):
    for line in fin[i]:
        xline= line.strip().split('\t')
        chrom= xline[0].split('_')[0]
        if chrom != prevChrom:
            if prevChrom:
                fout.close()
            ## Create output dir and open output file
            chromDir= os.path.join(outdir, chrom)
            if not os.path.exists(chromDir):
                 os.makedirs(chromDir)
            outfname= os.path.join(chromDir, os.path.split(infiles[i])[1])            
            fout= open(outfname, 'w')
            fout.write(headers[i] + '\n')
            fout.write(line)
            print('File "%s" - Writing to "%s"' %(infiles[i], outfname))            
            if not os.path.isdir(chromDir):
                try:
                    os.makedirs(chromDir)
                    ndirs += 1
                except OSError:
                    sys.exit('Cannot create dir %s' %(chromDir))
            prevChrom= chrom
        else:
            fout.write(line)
        if ndirs > 50:
            ## This is too make sure you don't flood the cwd with dirs and files!
            sys.exit('Too many dirs!')
    if ndirs > 50:
        sys.exit('Too many dirs!')
fout.close()
for f in fin:
    f.close()
