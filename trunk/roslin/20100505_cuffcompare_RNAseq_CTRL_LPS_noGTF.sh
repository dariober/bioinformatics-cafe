#! /bin/bash

# -----------------------------------------------------------------------------
#  Use cuffcompare to determine differenes between
#  RNAseq reads assembled by cufflinks w/o reference
#  GTF and ensembl transcriptome ss9.56
# -----------------------------------------------------------------------------

# ------------------------------[ User's input ]-------------------------------

# Output dir
cuffcomp_output='/exports/work/vet_roslin_nextgen/dario/cufflink/output/20100505_RNAseq_cuffcompare'

# Input GTF file(s) from cufflinks
gtf_1='/exports/work/vet_roslin_nextgen/dario/cufflink/output/20100317_RNAseq_LPS_noGTF/transcripts.gtf'
gtf_2='/exports/work/vet_roslin_nextgen/dario/cufflink/output/20100317_RNAseq_CTRL_noGTF/transcripts.gtf'

# Output prefix (-o)
outprefix='RNAseq'

# Optional reference GTF (-r)
ref_gtf='/exports/work/vet_roslin_nextgen/dario/GTF/Sus_scrofa.Sscrofa9.56.gtf'

# Note: output files .tmap and <transcript>.ref and <transcript>.refmap will be
# in the same dir as the input <transcript>.gtf

# -----------------------------[ Compile command ]-----------------------------

cd $cuffcomp_output

## /exports/work/vet_roslin_nextgen/dario/cufflink/current/cuffcompare -r $ref_gtf $ref_gtf $ref_gtf

/exports/work/vet_roslin_nextgen/dario/cufflink/current/cuffcompare -o $outprefix -r $ref_gtf -R $gtf_1 $gtf_2

