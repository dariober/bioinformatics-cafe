#!/bin/bash

# -----------------------------------------------------------------------------
# Prepare annotation database for pig
# Note: It may require nn GB of memory
# -----------------------------------------------------------------------------

# qsub -P vet_roslin -l h_rt=01:00:00 -pe memory 2 prepare_database.sh

source ~/.bash_profile

cd /exports/work/vet_roslin_nextgen/dario/annovar/database

# -----------------------------------------------------------------------------
# Un-comment as necessary:
#

## Sscrofa9/susScr2: ENSEMBL genes:
annotate_variation.pl -downdb -buildver susScr2 ensgene sscrofa2db
annotate_variation.pl --buildver susScr2 --downdb seq sscrofa2db/susScr2_seq
retrieve_seq_from_fasta.pl sscrofa2db/susScr2_ensGene.txt -seqdir sscrofa2db/susScr2_seq -format ensGene -outfile sscrofa2db/susScr2_ensGeneMrna.fa
