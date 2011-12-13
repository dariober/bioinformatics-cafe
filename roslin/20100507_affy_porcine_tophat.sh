#!/bin/bash

# 
# Align probes from affymetrix porcine chip na30 against pig genome using tophat
# 

#----------------------[ Add bowtie and tophat to PATH ]-----------------------
# Uncomment as necessary
PATH=$PATH:/exports/work/vet_roslin_nextgen/dario/bowtie/current/
PATH=$PATH:/exports/work/vet_roslin_nextgen/dario/tophat/current/bin/

tophat \
    --output-dir /exports/work/vet_roslin_nextgen/dario/tophat/output/20100507_affy_porcine \
    --GFF /exports/work/vet_roslin_nextgen/dario/GFF/Sus_scrofa.Sscrofa9.56.gff3 \
    /exports/work/vet_roslin_nextgen/dario/bowtie/indexes/Sscrofa9.56_plus_Human_ribosomal_DNA_complete_repeating_unit \
    /exports/work/vet_roslin_nextgen/dario/fasta/affymetrix_porcine_na30/Porcine_na30.fa
  
exit

#  Default values
# --min-anchor-length  8
# --segment-length 25
# --max-multihits 40
# --splice-mismatches 0

  