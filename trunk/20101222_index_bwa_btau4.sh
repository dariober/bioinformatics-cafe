#!/bin/bash

#
#  Download B.tau genome 4.0 and build index for BWA
#
# Using:
# qsub -P vet_roslin -pe memory 2 -l h_rt=05:59:00 20101222_index_bwa_btau4.sh

## Add bwa to path
PATH=$PATH:/exports/work/vet_roslin_nextgen/dario/bwa/bwa-0.5.8c

## FASTA and index files here:
cd /exports/work/vet_roslin_nextgen/dario/bwa/indexes

# ---------------[ Download and unzip genome from ensembl ]--------------------

# wget ftp://ftp.ensembl.org/pub/current/fasta/bos_taurus/dna/Bos_taurus.Btau_4.0.60.dna.toplevel.fa.gz
# gunzip Bos_taurus.Btau_4.0.60.dna.toplevel.fa.gz
# -----------------------------------------------------------------------------

## Now build index:
bwa index -a bwtsw /exports/work/vet_roslin_nextgen/dario/bwa/indexes/Bos_taurus.Btau_4.0.60.dna.toplevel.fa