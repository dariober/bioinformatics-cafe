#!/bin/bash

cd /exports/work/vet_roslin_nextgen/dario/bwa/output/20101218_cage_050810_increment
PATH=$PATH:/exports/work/vet_roslin_nextgen/dario/samtools/samtools-0.1.10

## Pileup cage tags mapping to the forward strand
samtools pileup \
    -c \
    -s \
    -f /exports/work/vet_roslin_nextgen/dario/bwa/indexes/Sus_scrofa.Sscrofa9.56.dna.toplevel.fa \
    cage_050810_bwa_20101218.clean.single.forward.bam > cage_050810_bwa_20101218.clean.single.forward.pileup

## Pileup cage tags mapping to the revesre strand
samtools pileup \
    -c \
    -s \
    -f /exports/work/vet_roslin_nextgen/dario/bwa/indexes/Sus_scrofa.Sscrofa9.56.dna.toplevel.fa \
    cage_050810_bwa_20101218.clean.single.reverse.bam > cage_050810_bwa_20101218.clean.single.reverse.pileup