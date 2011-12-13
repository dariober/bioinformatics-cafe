#!/bin/bash

# -----------------------------------------------------------------------------
# Use cuffcompare to compare libraries RNAseq_LPS vs RNAseq_CTRL vs pig 
# GTF annotation. 
# Note: Version 9.57 of gtf was released but not used (use 9.56 to comply with
# tophat and cufflinks analyses)
# -----------------------------------------------------------------------------

## Where output will be:
mkdir /exports/work/vet_roslin_nextgen/dario/cufflink/output/20100316_cuffcompare_RNAseq/

cd /exports/work/vet_roslin_nextgen/dario/cufflink/output/20100316_cuffcompare_RNAseq/

/exports/work/vet_roslin_nextgen/dario/cufflink/current/cuffcompare \
  /exports/work/vet_roslin_nextgen/dario/cufflink/output/20100301_RNAseq_CTRL/transcripts.gtf \
  /exports/work/vet_roslin_nextgen/dario/cufflink/output/20100301_RNAseq_LPS/transcripts.gtf \
  -o 20100316_cuffcompare_RNAseq_summary.txt \
  -r /exports/work/vet_roslin_nextgen/dario/GTF/Sus_scrofa.Sscrofa9.56.gtf \
  -R \
  -V 