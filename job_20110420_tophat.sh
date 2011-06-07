#!/bin/bash

#
#  Submit with 
#  qsub -P vet_roslin -l h_rt=23:00:00 -pe memory 4 job_20110420_tophat.sh
#

PATH=$PATH:/exports/work/vet_roslin_nextgen/dario/tophat/tophat-1.2.0.Linux_x86_64
PATH=$PATH:/exports/work/vet_roslin_nextgen/dario/samtools/samtools-0.1.15
PATH=$PATH:/exports/work/vet_roslin_nextgen/dario/bowtie/bowtie-0.12.7

## Working dir (output will go here)
cd /exports/work/vet_roslin_nextgen/dario/tophat/output/20110420_rnaseq_am_bmdm/bmdm/ctrl

## Path to where reads are
read_path="/exports/work/vet_roslin_nextgen/dario/fastq/20100614_RNAseq_pig_BMDM/478"

tophat \
  --output-dir ./ \
  -r 100 \
  --mate-std-dev 100 \
  --num-threads 4 \
  -G /exports/work/vet_roslin_nextgen/dario/ensembl/release-60/gtf/sus_scrofa/Sus_scrofa.Sscrofa9.60.gtf \
  --allow-indels \
  /exports/work/vet_roslin_nextgen/dario/bowtie/indexes/Sus_scrofa.Sscrofa9.56.dna.toplevel \
  $read_path/BMDM_CTRL_14JUN10_300BP_s_7_1_sequence.sanger.fq.trimmed $read_path/BMDM_CTRL_14JUN10_300BP_s_7_2_sequence.sanger.fq.trimmed

samtools index accepted_hits.bam
samtools flagstat accepted_hits.bam

exit