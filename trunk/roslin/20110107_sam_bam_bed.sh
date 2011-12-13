#!/bin/bash

#
# Convert SAM file to BAM and this latter to bed.
#


## Bedtools and SAMtools
PATH=$PATH:/exports/work/vet_roslin_nextgen/dario/BEDTools/BEDTools-Version-2.10.1/bin
PATH=$PATH:/exports/work/vet_roslin_nextgen/dario/samtools/samtools-0.1.10

cd /exports/work/vet_roslin_nextgen/dario/bwa/output/20110106_fantom5_macro_monocyte_vs_ss9

## If necessary prepare the .fai header with 
## samtools faidx reference.fasta

## Convert SAM to BAM using the appropriate .fai header:
echo 'Converting to BAM ...'
samtools view -b -S \
    -o tmp_unsorted.bam \
    -t /exports/work/vet_roslin_nextgen/dario/bwa/indexes/Sus_scrofa.Sscrofa9.56.dna.toplevel.fa.fai \
    mphage_mcyte_f5_hg19_21bp.nordna.single.sam.gz

## Sort BAM:
echo 'Sorting BAM ...'
samtools sort tmp_unsorted.bam mphage_mcyte_f5_hg19_21bp.nordna.single
rm tmp_unsorted.bam

## Index as well
echo 'Indexing BAM...'
samtools index mphage_mcyte_f5_hg19_21bp.nordna.single.bam

## Convert BAM to BED:
echo 'Converting BAM to BED...'
bamToBed -i mphage_mcyte_f5_hg19_21bp.nordna.single.bam > mphage_mcyte_f5_hg19_21bp.nordna.single.bed