#!/bin/bash

# -----------------------------------------------------------------------------
# Map RNAseq libraries against pig genome + human rDNA.
# This output used for allele specific expression.
# Report one match per read (-m 1), this match being from the best stratum (--best --strata). 
# Suppress reads that align in more than 1 position.
# Reads filtered to have at least 15 bp with no more than 1 base with qs<20.
# Using bowtie-0.12.4.
# -----------------------------------------------------------------------------

## Where output will be:
/exports/work/vet_roslin_nextgen/dario/bowtie/current/bowtie \
  /exports/work/vet_roslin_nextgen/dario/bowtie/indexes/Sscrofa9.56_plus_Human_ribosomal_DNA_complete_repeating_unit \
  -1 /exports/work/vet_roslin_nextgen/dario/fastq/20091217_RNAseq_AM_LPS_CTRL/filtered_15bp_phred20/s_7_1_sequence.txt.filtered \
  -2 /exports/work/vet_roslin_nextgen/dario/fastq/20091217_RNAseq_AM_LPS_CTRL/filtered_15bp_phred20/s_7_2_sequence.txt.filtered \
  -q \
  --solexa1.3-quals \
  -a \
  -m 1 \
  -v 2 \
  --best \
  --strata \
  --sam \
  /exports/work/vet_roslin_nextgen/dario/bowtie/output/20100406_RNAseq_LPS_sscrofa9.56/20100406_RNAseq_LPS_sscrofa9.56.sam
