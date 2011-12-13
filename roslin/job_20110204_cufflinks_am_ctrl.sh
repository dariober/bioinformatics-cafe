#!/bin/bash

# -----------------------------------------------------------------------------
# Execute cufflinks using GTF reference (no denovo assembly)
# -----------------------------------------------------------------------------


PATH=$PATH:/exports/work/vet_roslin_nextgen/dario/cufflinks/cufflinks-0.9.3.Linux_x86_64

cd /exports/work/vet_roslin_nextgen/dario/tophat/output/20110202_rnaseq_am_ctrl

cufflinks   accepted_hits.bam \
  -G /exports/work/vet_roslin_nextgen/dario/ensembl/release-60/gtf/sus_scrofa/Sus_scrofa.Sscrofa9.60.gtf \
  -r /exports/work/vet_roslin_nextgen/dario/bowtie/indexes/Sus_scrofa.Sscrofa9.56.dna.toplevel.fa \
  -N \ 

mv transcripts.expr transcripts_gtf.expr
mv transcripts.gtf transcripts_gtf.gtf
