#!/bin/bash

# -----------------------------------------------------------------------------
# Align RNAseq tags against pig transcript ENSSSCT00000005341.
# This transcript is supposed to be the pig homolog of human CHRNA7 (see lab 
# book 6/10/2010 and later)
# -----------------------------------------------------------------------------

cd /exports/work/vet_roslin_nextgen/dario/bowtie/output/20101008_chrna7

/exports/work/vet_roslin_nextgen/dario/bowtie/current/bowtie \
    /exports/work/vet_roslin_nextgen/dario/bowtie/indexes/ENSSSCT00000005341 \
    -1 /exports/work/vet_roslin_nextgen/dario/fastq/20091217_RNAseq_AM_LPS_CTRL/s_7_1_sequence.txt \
    -2 /exports/work/vet_roslin_nextgen/dario/fastq/20091217_RNAseq_AM_LPS_CTRL/s_7_2_sequence.txt \
    -q \
    -a \
    -m 10 \
    -v 2 \
    --best \
    --strata \
    --sam \
    /exports/work/vet_roslin_nextgen/dario/bowtie/output/20101008_chrna7/ENSSSCT00000005341_lps.sam