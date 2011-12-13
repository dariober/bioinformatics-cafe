#!/bin/bash

#
#  Make index for bowtie of the S. scrofa genome, version 9.56
#


cd /exports/work/vet_roslin_nextgen/dario/bowtie/indexes
/exports/work/vet_roslin_nextgen/dario/bowtie/current/bowtie-build \
    /exports/work/vet_roslin_nextgen/dario/genomes/sscrofa9.56/Sus_scrofa.Sscrofa9.56.dna.toplevel.fa \
    Sus_scrofa.Sscrofa9.56.dna.toplevel 

exit


## Test: make index of chr 11

cd /exports/work/vet_roslin_nextgen/dario/bowtie/indexes
/exports/work/vet_roslin_nextgen/dario/bowtie/current/bowtie-build \
    /exports/work/vet_roslin_nextgen/dario/genomes/sscrofa9.56/chromosomes/Sus_scrofa.Sscrofa9.56.dna_rm.chromosome.18.fa \
    Sus_scrofa.Sscrofa9.56.dna_rm.chromosome.18 

