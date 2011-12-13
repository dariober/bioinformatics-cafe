#!/bin/bash

#
#  Align unique CAGE tags from humans against sscrofa9.56 (with human rDNA).
#  CAGE tags obtained from FANTOM4 (see LabBook on 16/08/2010)
#
#  Report alignments where the best match (--best) in the best stratum (--strata) is unique (-m 1).
#

cd /exports/work/vet_roslin_nextgen/dario/bowtie/output/20100816_human_CAGE_vs_sscrofa9

/exports/work/vet_roslin_nextgen/dario/bowtie/current/bowtie \
    /exports/work/vet_roslin_nextgen/dario/bowtie/indexes/Sscrofa9.56_plus_Human_ribosomal_DNA_complete_repeating_unit \
    -f /exports/work/vet_roslin_nextgen/dario/fasta/human_cage_tags/human_cage_tags_unique.fa \
    -a \
    -m 1 \
    -v 2 \
    --best \
    --strata \
    --sam \
    -p 3 \
    /exports/work/vet_roslin_nextgen/dario/bowtie/output/20100816_human_CAGE_vs_sscrofa9/20100816_human_cage_vs_sscrofa9.sam