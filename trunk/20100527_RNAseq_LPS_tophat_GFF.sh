#!/bin/bash

# Use TopHat 1.0.13 to align reads from RNAseq_CTRL (see 24/5/2010 for filtering) against sscrofa 9.56
# Supply a GFF annotation file.
# Note the use of --best --strata --seedmms (not allowed in original script /bin/tophat)

#----------------------[ Add bowtie and tophat to PATH ]-----------------------
# Uncomment as necessary
PATH=$PATH:/exports/work/vet_roslin_nextgen/dario/bowtie/current/
PATH=$PATH:/exports/work/vet_roslin_nextgen/dario/tophat/current/bin/

tophat \
  --output-dir /exports/work/vet_roslin_nextgen/dario/tophat/output/20100527_RNAseq_LPS/ \
  --mate-inner-dist 130 \
  --mate-std-dev 30 \
  --solexa1.3-quals \
  --GFF /exports/work/vet_roslin_nextgen/dario/GFF/Sus_scrofa.Sscrofa9.56.gff3 \
  --best \
  --strata \
  --seedmms 3 \
  /exports/work/vet_roslin_nextgen/dario/bowtie/indexes/Sus_scrofa.Sscrofa9.56.dna.toplevel \
  /exports/work/vet_roslin_nextgen/dario/fastq/20091217_RNAseq_AM_LPS_CTRL/20100527_filter_LPS/s_7_1_sequence.txt.avgq \
  /exports/work/vet_roslin_nextgen/dario/fastq/20091217_RNAseq_AM_LPS_CTRL/20100527_filter_LPS/s_7_2_sequence.txt.avgq

#  Default values
# --min-anchor-length  8
# --segment-length 25
# --max-multihits 40
# --splice-mismatches 0
  