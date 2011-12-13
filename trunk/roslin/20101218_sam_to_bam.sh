#!/bin/bash

cd /exports/work/vet_roslin_nextgen/dario/bwa/output/20101218_cage_050810_increment
PATH=$PATH:/exports/work/vet_roslin_nextgen/dario/samtools/samtools-0.1.10


## If necessary prepare the .fai header with 
## samtools faidx reference.fasta

## Convert SAM to BAM using the appropriate .fai header:
echo 'Converting to BAM ...'
samtools view -b -S \
    -o tmp_unsorted.bam \
    -t /exports/work/vet_roslin_nextgen/dario/bwa/indexes/Sus_scrofa.Sscrofa9.56.dna.toplevel.fa.fai \
    cage_050810_bwa_20101218.clean.single.sam

## Sort BAM:
echo 'Sorting BAM ...'
samtools sort tmp_unsorted.bam cage_050810_bwa_20101218.clean.single
rm tmp_unsorted.bam

## Extract forward and reverse mapping reads:
## Forward:
echo 'Extracting forward reads...'
samtools view -b -F 16 -o cage_050810_bwa_20101218.clean.single.forward.bam cage_050810_bwa_20101218.clean.single.bam
samtools index cage_050810_bwa_20101218.clean.single.forward.bam

# Optional re-convert to SAM
# samtools view -o cage_050810_bwa_20101218.single.forward.sam cage_050810_bwa_20101218.clean.single.forward.bam

## Reverse:
echo 'Extracting reverse reads...'
samtools view -b -f 16 -o cage_050810_bwa_20101218.clean.single.reverse.bam cage_050810_bwa_20101218.clean.single.bam
samtools index cage_050810_bwa_20101216.clean.single.reverse.bam

# Optional re-convert to SAM
# samtools view -o cage_050810_bwa_20101218.clean.single.reverse.sam cage_050810_bwa_20101218.clean.single.reverse.bam
