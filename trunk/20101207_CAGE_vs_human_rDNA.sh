#!/bin/bash

# -----------------------------------------------------------------------------
# Map CAGE library on labbbok 06/12/2010 against pig genome + human rDNA.
# Report up to 10 matches per read (-m 10) from the error best stratum (--best --strata).
# Suppress reads that align in more than 10 positions.
# Reads trimmed 9 bases from 5' end (--trim5) to remove AGACAGCAG
# -----------------------------------------------------------------------------

cd /exports/work/vet_roslin_nextgen/dario/bowtie/output/20101207_CAGE_vs_sscrofa9_human_rDNA


## Where output will be:
/exports/work/vet_roslin_nextgen/dario/bowtie/bowtie-0.12.7/bowtie \
  /exports/work/vet_roslin_nextgen/dario/bowtie/indexes/Sscrofa9.56_plus_Human_ribosomal_DNA_complete_repeating_unit \
  /exports/work/vet_roslin_nextgen/dario/fastq/20100805_CAGE_AM/Beraldi_CAGE_50810_101202_EBRI093151_0060_s_8_sequence.txt  -q \
  --solexa1.3-quals \
  --trim5 9 \
  -a \
  -m 10 \
  --maqerr 70 \
  --nomaqround \
  --best \
  --strata \
  --sam \
  /exports/work/vet_roslin_nextgen/dario/bowtie/output/20101207_CAGE_vs_sscrofa9_human_rDNA/20101207_CAGE_vs_sscrofa9_human_rDNA.sam
