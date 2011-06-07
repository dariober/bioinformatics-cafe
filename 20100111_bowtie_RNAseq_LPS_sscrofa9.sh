## Align pair-end reads from RNAseq library, LPS treated, received 17/12/2009
## Reads filtered see 17/12/2009
## Alignement against S scrofa 9.56 with human rDNA included
## Allow 2 mismatches
##

#!/usr/bin/bash

/exports/work/vet_roslin_nextgen/dario/bowtie/bowtie-0.12.1/bowtie \
  Sscrofa9.56_plus_Human_ribosomal_DNA_complete_repeating_unit \
  -1 /exports/work/vet_roslin_nextgen/dario/fastq/20091217_RNAseq_AM_LPS_CTRL/filtered/s_7_1_sequence.txt.filtered \
  -2 /exports/work/vet_roslin_nextgen/dario/fastq/20091217_RNAseq_AM_LPS_CTRL/filtered/s_7_2_sequence.txt.filtered \
  -q \
  --solexa1.3-quals \
  -M 10 \
  -v 2 \
  --sam \
  /exports/work/vet_roslin_nextgen/dario/bowtie/output/bwt_RNAseq_AM_LPS.map
  