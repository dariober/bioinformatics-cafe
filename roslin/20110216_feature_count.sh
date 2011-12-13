#!/bin/bash

#
# Count reads mapping to a GTF feature.
# Input: bam file (typically from tophat)
# Output: Output of HTSeq
#

PATH=$PATH:/exports/work/vet_roslin_nextgen/dario/samtools/samtools-0.1.10

module add python/2.6.3

## Working dir
cd /exports/work/vet_roslin_nextgen/dario/tophat/output/20110202_rnaseq_am_ctrl/

## Convert bam to sam:
samtools view accepted_hits.bam > accepted_hits.sam

## Sort sam file by read name (actually by entire line)
python /exports/work/vet_roslin_nextgen/dario/bin/sort_file-2.4.py accepted_hits.sam accepted_hits.qnamesorted.sam

## Count reads per feature
python -m HTSeq.scripts.count -m union -s no -t exon -i gene_id \
    accepted_hits.qnamesorted.sam \
    /exports/work/vet_roslin_nextgen/dario/ensembl/release-60/gtf/sus_scrofa/Sus_scrofa.Sscrofa9.60.gtf > feature_count.htseq
    