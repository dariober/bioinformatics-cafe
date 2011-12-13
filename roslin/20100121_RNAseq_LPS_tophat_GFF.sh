#!/bin/bash

# Use TopHat 1.0.12 to align reads from RNAseq_LPS (only reads 35bp long
# and no more than two bases with phred < 20) against sscrofa 9.56 + human rDNA
# Supply a GFF annotation file.


#----------------------[ Add bowtie and tophat to PATH ]-----------------------
# Uncomment as necessary
PATH=$PATH:/exports/work/vet_roslin_nextgen/dario/bowtie/current/
PATH=$PATH:/exports/work/vet_roslin_nextgen/dario/tophat/current/bin/

tophat \
  --output-dir /exports/work/vet_roslin_nextgen/dario/tophat/output/20100121_RNAseq_LPS/ \
  --mate-inner-dist 130 \
  --mate-std-dev 30 \
  --solexa1.3-quals \
  --GFF /exports/work/vet_roslin_nextgen/dario/GFF/Sus_scrofa.Sscrofa9.56.gff3.ensnames \
  /exports/work/vet_roslin_nextgen/dario/bowtie/current/indexes/Sscrofa9.56_plus_Human_ribosomal_DNA_complete_repeating_unit \
  /exports/work/vet_roslin_nextgen/dario/fastq/20091217_RNAseq_AM_LPS_CTRL/filtered_35bp_phred20/RNAseq_LPS_1_35bp.fastq \
  /exports/work/vet_roslin_nextgen/dario/fastq/20091217_RNAseq_AM_LPS_CTRL/filtered_35bp_phred20/RNAseq_LPS_2_35bp.fastq \

#  Default values
# --min-anchor-length  8
# --segment-length 25
# --max-multihits 40
# --splice-mismatches 0

  