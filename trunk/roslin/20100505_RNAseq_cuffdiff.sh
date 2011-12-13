#!/bin/bash

# -----------------------------------------------------------------------------
# 
#
# -----------------------------------------------------------------------------

## Where output will be:
outdir='/exports/work/vet_roslin_nextgen/dario/cufflink/output/20100505_RNAseq_cuffdiff/'

gtf='/exports/work/vet_roslin_nextgen/dario/cufflink/output/20100505_RNAseq_cuffcompare/RNAseq.combined.gtf'
sam1='/exports/work/vet_roslin_nextgen/dario/tophat/output/20100122_RNAseq_CTRL/accepted_hits.sam'
sam2='/exports/work/vet_roslin_nextgen/dario/tophat/output/20100122_RNAseq_LPS/accepted_hits.sam'

# mkdir $outdir
cd $outdir

/exports/work/vet_roslin_nextgen/dario/cufflink/current/cuffdiff $gtf $sam1 $sam2 --inner-dist-mean 130 -s 30