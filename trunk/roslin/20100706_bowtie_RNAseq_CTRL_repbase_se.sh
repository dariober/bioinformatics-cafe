#!/bin/bash

#
#  Align RNAseq_CTRL and RNAseq_LPS libs 
#  - Note the use of the index with reference names renamed
#  - Note that pair end reads are run as individual single end reads with the -c option
#

cd /exports/work/vet_roslin_nextgen/dario/bowtie/output/

/exports/work/vet_roslin_nextgen/dario/bowtie/current/bowtie \
  /exports/work/vet_roslin_nextgen/dario/bowtie/indexes/pig_repbase_20090604_renamed \
  /exports/work/vet_roslin_nextgen/dario/fastq/20091217_RNAseq_AM_LPS_CTRL/s_8_1_sequence.txt,\
/exports/work/vet_roslin_nextgen/dario/fastq/20091217_RNAseq_AM_LPS_CTRL/s_8_2_sequence.txt \
  -q \
  --solexa1.3-quals \
  -a \
  --seedmms 3 \
  --best \
  --strata \
  --sam \
  --threads 2 \
  /exports/work/vet_roslin_nextgen/dario/bowtie/output/20100706_RNAseq_repbase_se/20100706_RNAseq_CTRL_repbase_se.sam