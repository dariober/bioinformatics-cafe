#!/bin/bash

# 
# Count reads in BAM file mapping to a gtf feature (gene id) using BEDTools/coverageBed.
# Remove multi-mappers first.
# 

## Note: Tophat doesn't really give MAPQ quality. Anything > 3 means single mapper.

source ~/.bash_profile
gtf='/exports/work/vet_roslin_nextgen/dario/ensembl/release-60/gtf/sus_scrofa/Sus_scrofa.Sscrofa9.60.gtf'

# -----------------------------------------------------------------------------
cd /exports/work/vet_roslin_nextgen/dario/tophat/output/20110202_rnaseq_am_ctrl
samtools view -q 15 -b accepted_hits.bam | coverageBed -abam stdin -b $gtf > accepted_hits.coverageBed

# -----------------------------------------------------------------------------
cd /exports/work/vet_roslin_nextgen/dario/tophat/output/20110202_rnaseq_am_lps
samtools view -q 15 -b accepted_hits.bam | coverageBed -abam stdin -b $gtf > accepted_hits.coverageBed

# -----------------------------------------------------------------------------
cd /exports/work/vet_roslin_nextgen/dario/tophat/output/20110202_rnaseq_bmdm_ctrl
samtools view -q 15 -b accepted_hits.bam | coverageBed -abam stdin -b $gtf > accepted_hits.coverageBed

# -----------------------------------------------------------------------------
cd /exports/work/vet_roslin_nextgen/dario/tophat/output/20110202_rnaseq_bmdm_lps
samtools view -q 15 -b accepted_hits.bam | coverageBed -abam stdin -b $gtf > accepted_hits.coverageBed