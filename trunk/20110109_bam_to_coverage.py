#!/usr/bin/python

""" Compute coverage of a the feature in GTF file for a list of BAM files.
    BAM files are filtered for a given MAPQ score:
.../samtools-0.1.10/samtools view -b -h -q 15 LN_173.fq.bam > tmp_filtered.bam
.../BEDTools-Version-2.10.1/bin/coverageBed -abam tmp_filtered.bam -b /.../bos_taurus/Bos_taurus.Btau_4.0.60.gtf> LN_173.cov
    SAM files have header. (if not use "samtools faidx my_reference.fasta")

"""

## Path to samtools / bedtools
beddir= '/exports/work/vet_roslin_nextgen/dario/BEDTools/BEDTools-Version-2.10.1/bin'
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

samtools= os.path.join(samdir, 'samtools view -b -h -q %d '%q )
bedtools= os.path.join(beddir, 'coverageBed ' )
for bam in bamfiles:
    " Filter for MAPQ "
    fcmd= samtools + bam  + ' > tmp_filtered.bam'
    print('Filtering with: ' + fcmd )
    os.system(fcmd)

    " Coverage "
    cov= re.sub('.fq.bam$', '.cov', bam) ## Output file
    covcmd= bedtools + '-abam tmp_filtered.bam -b %s> '%gtf + cov 
    print('Executing: ' + covcmd)
    os.system(covcmd)
os.remove('tmp_filtered.bam')

