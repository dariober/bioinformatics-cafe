#!/bin/bash

# -----------------------------------------------------------------------------
# Isoform and transcript quantitation analysis of RNAseq_CTRL using cufflink 
# Input alignment file produced by tophat 
# -----------------------------------------------------------------------------

## Where output will be:
mkdir /exports/work/vet_roslin_nextgen/dario/cufflink/output/20100317_RNAseq_CTRL_noGTF/

cd /exports/work/vet_roslin_nextgen/dario/cufflink/output/20100317_RNAseq_CTRL_noGTF/

/exports/work/vet_roslin_nextgen/dario/cufflink/current/cufflinks \
  /exports/work/vet_roslin_nextgen/dario/tophat/output/20100122_RNAseq_CTRL/accepted_hits.sam \
  --inner-dist-mean 130  \
  -s 30 \
  -L CTRL \  
  