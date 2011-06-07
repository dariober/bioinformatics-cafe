#!/bin/bash

cd /exports/work/vet_roslin_nextgen/dario/bwa/output/20101216_cage_050810_iterate
PATH=$PATH:/exports/work/vet_roslin_nextgen/dario/samtools/samtools-0.1.10

## Pileup cage tags mapping to the forward strand
samtools pileup \
    -c \
    -s \
    -f /exports/work/vet_roslin_nextgen/dario/bwa/indexes/Sscrofa9.56_plus_Human_ribosomal_DNA_complete_repeating_unit.fa \
    cage_050810_bwa_20101216.clean.single.forward.bam > cage_050810_bwa_20101216.clean.single.forward.pileup

## Pileup cage tags mapping to the revesre strand
samtools pileup \
    -c \
    -s \
    -f /exports/work/vet_roslin_nextgen/dario/bwa/indexes/Sscrofa9.56_plus_Human_ribosomal_DNA_complete_repeating_unit.fa \
    cage_050810_bwa_20101216.clean.single.reverse.bam > cage_050810_bwa_20101216.clean.single.reverse.pileup