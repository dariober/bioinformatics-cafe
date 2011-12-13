#!/usr/bin/python

""" Compute coverage of a the feature in GTF file for a list of BAM files.
    Using python -m HTSeq.scripts.count
    Alignments are filtered for a given MAPQ score:
    " samtools view -q 15 LN_109.fq.bam > tmp.sam "
    " python -m HTSeq.scripts.count -s -m union -t exon -i gene_id tmp.sam Bos_taurus.Btau_4.0.60.gtf > LN_109.cov"

"""

## Issue module add python/2.6.3 before starting to use python 2.6.3

## Path to samtools / bedtools
samdir= '/exports/work/vet_roslin_nextgen/dario/samtools/samtools-0.1.10/'

## List of bam files to compute coverage
bamfiles= ['LN_109.fq.bam', 'LN_114.fq.bam', 'LN_124.fq.bam', 'LN_130.fq.bam', 'LN_146.fq.bam', 'LN_173.fq.bam', 'LN_183.fq.bam', 'LN_20B.fq.bam', 'LN_21.fq.bam', 'LN_38.fq.bam', 'LN_47.fq.bam', 'LN_50.fq.bam', 'LN_57.fq.bam', 'LN_58.fq.bam', 'LN_92.fq.bam']

## GTF file to use to compute coverage
gtf= '/exports/work/vet_roslin_nextgen/dario/ensembl/release-60/gtf/bos_taurus/Bos_taurus.Btau_4.0.60.gtf'

## Mapping quality threshols (samtools -q):
q= 15

## Output is the same as input with .bam replaced by .cov

# -----------------------------------------------------------------------------

import os
import re

samtools= os.path.join(samdir, 'samtools view -h -q %d '%q )
 
for bam in bamfiles:
    " Filter for MAPQ "
    fcmd= samtools + bam  + ' > tmp_filtered.sam'
    print('Filtering with: ' + fcmd )
    os.system(fcmd)

    " Coverage "
    cov= re.sub('.fq.bam$', '.cov', bam) ## Output file
    htseq= 'python -m HTSeq.scripts.count -q -s yes -m union -t exon -i gene_id tmp_filtered.sam ' + gtf + ' > ' + cov
    print('Executing: ' + htseq)
    os.system(htseq)
os.remove('tmp_filtered.sam')

