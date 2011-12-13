#!/bin/bash

# Submit with
# qsub -P vet_roslin -pe memory 4 -l h_rt=05:59:00 job_20110121_test_aln.sh

## Add bwa to PATH
PATH=$PATH:/exports/work/vet_roslin_nextgen/dario/bwa/bwa-0.5.8c

## Make a directory for this alignment
## mkdir /exports/work/vet_roslin_nextgen/dario/bwa/output/20110121_mouse_chipseq_mb
cd /exports/work/vet_roslin_nextgen/dario/bwa/output/20110121_mouse_chipseq_mb

bwa aln -t 4 /exports/work/vet_roslin_nextgen/dario/bwa/indexes/Mus_musculus.NCBIM37.60.dna.toplevel.fa /exports/work/vet_roslin_nextgen/dario/Tritume/GSM611116_2060_s_4_export.fq > GSM611116_2060_s_4_export.bwa

bwa samse /exports/work/vet_roslin_nextgen/dario/bwa/indexes/Mus_musculus.NCBIM37.60.dna.toplevel.fa GSM611116_2060_s_4_export.bwa /exports/work/vet_roslin_nextgen/dario/Tritume/GSM611116_2060_s_4_export.fq > GSM611116_2060_s_4_export.sam
