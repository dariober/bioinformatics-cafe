#!/bin/bash

# -----------------------------------------------------------------------------
# Isoform and transcript quantitation analysis of RNAseq_CTRL using cufflink 
# Input alignment file produced by tophat 
# -----------------------------------------------------------------------------

## Where output will be:
out='/exports/work/vet_roslin_nextgen/dario/cufflink/output/20100602_RNAseq_CTRL2_denovo/'

mkdir $out
cd $out

/exports/work/vet_roslin_nextgen/dario/cufflink/current/cufflinks \
  /exports/work/vet_roslin_nextgen/dario/tophat/output/20100528_RNAseq_CTRL/accepted_hits.sam \
  --inner-dist-mean 130  \
  -s 30 \
  -L CTRL2_denovo  
  