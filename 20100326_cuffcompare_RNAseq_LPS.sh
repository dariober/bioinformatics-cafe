#! /bin/bash

# -----------------------------------------------------------------------------
#  Use cuffcompare to determine differenes between
#  RNAseq LPS assembled by cufflinks w/o reference GTF 
#  Ensembl transcriptome ss9.56
# -----------------------------------------------------------------------------

# ------------------------------[ User's input ]-------------------------------

# Output dir
cuffcomp_output='/exports/work/vet_roslin_nextgen/dario/cufflink/output/20100326_cuffcompare_RNAseq_LPS_noGTF'

# Input GTF file(s) from cufflinks
gtf_1='/exports/work/vet_roslin_nextgen/dario/cufflink/output/20100317_RNAseq_LPS/transcripts.gtf'
gtf_2='/exports/work/vet_roslin_nextgen/dario/cufflink/output/20100317_RNAseq_CTRL/transcripts.gtf'

# Output prefix (-o)
outprefix='RNAseq'

# Optional reference GTF (-r)
ref_gtf='/exports/work/vet_roslin_nextgen/dario/GTF/test_Sus_scrofa.Sscrofa9.56.gtf'

# Note: output files .tmap and <transcript>.ref and <transcript>.refmap will be
# in the same dir as the input <transcript>.gtf

# -----------------------------[ Compile command ]-----------------------------

mkdir $cuffcomp_output
cd $cuffcomp_output

/exports/work/vet_roslin_nextgen/dario/cufflink/current/cuffcompare -o $outprefix -r $ref_gtf -R $gtf_1

