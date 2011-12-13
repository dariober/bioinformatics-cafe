## Align CAGE inserts from libraries #26 and #27 cloned and sequenced.
## See 23/02/2010 and earlier.
## Alignement against S scrofa 9.56 with human rDNA included
## Allow 2 mismatches
## Report up to 10 hits per read. For reads with more than 10 hits, report one hit at random.

#!/usr/bin/bash

/exports/work/vet_roslin_nextgen/dario/bowtie/current/bowtie \
  Sscrofa9.56_plus_Human_ribosomal_DNA_complete_repeating_unit \
  /exports/work/vet_roslin_nextgen/dario/bowtie/current/reads/cage_inserts_20100223.fa \
  -f \
  -a \
  --best \
  --strata \
  -m 10 \
  -v 3 \
  --sam \
##  /exports/work/vet_roslin_nextgen/dario/bowtie/output/bwt_20100223_cloned_inserts.sam
  