#!/bin/bash

#
# Convert SAM file to BAM and this latter to bed.
# Here we convert CAGE tags 050810 mapped with bwa (incerment strategy), filtered 
# for rDNA contamination, single mappers only.
# The BED file of aligned CAGE tags gives the transcription start sites ('feature start'
# from (+)tags and 'feature end' for (-)tags)


## Bedtools and SAMtools
PATH=$PATH:/exports/work/vet_roslin_nextgen/dario/BEDTools/BEDTools-Version-2.10.1/bin
PATH=$PATH:/exports/work/vet_roslin_nextgen/dario/samtools/samtools-0.1.10

cd /exports/work/vet_roslin_nextgen/dario/bwa/output/20101218_cage_050810_increment

## If necessary prepare the .fai header with 
## samtools faidx reference.fasta

## Convert SAM to BAM using the appropriate .fai header:
#  echo 'Converting to BAM ...'
#  samtools view -b -S \
#    -o cage_050810_bwa_20101218.clean.single.bam \
#    -t /exports/work/vet_roslin_nextgen/dario/bwa/indexes/Sus_scrofa.Sscrofa9.56.dna.toplevel.fa.fai \
#    cage_050810_bwa_20101218.clean.single.sam

## Convert BAM to BED:
echo 'Converting BAM to BED...'
bamToBed -i cage_050810_bwa_20101218.clean.single.bam > cage_050810_bwa_20101218.nordna.single.bed