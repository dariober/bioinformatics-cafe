#!/bin/bash

#
#  Submit with 
#  qsub -P vet_roslin -l h_rt=23:00:00 -pe memory 4 job_20110202_rnaseq_am_ctrl.sh
#


PATH=$PATH:/exports/work/vet_roslin_nextgen/dario/tophat/tophat-1.2.0.Linux_x86_64
PATH=$PATH:/exports/work/vet_roslin_nextgen/dario/samtools/samtools-0.1.10
PATH=$PATH:/exports/work/vet_roslin_nextgen/dario/bowtie/bowtie-0.12.7

## Working dir (output will go here)
cd /exports/work/vet_roslin_nextgen/dario/tophat/output/20110202_rnaseq_am_ctrl

## Path to where reads are
read_path="/exports/work/vet_roslin_nextgen/dario/fastq/20091217_RNAseq_AM_LPS_CTRL"

tophat \
  --output-dir ./ \
  -r 100 \
  --mate-std-dev 100 \
  --solexa1.3-quals \
  --num-threads 4 \
  -G /exports/work/vet_roslin_nextgen/dario/ensembl/release-60/gtf/sus_scrofa/Sus_scrofa.Sscrofa9.60.gtf \
  --allow-indels \
  /exports/work/vet_roslin_nextgen/dario/bowtie/indexes/Sus_scrofa.Sscrofa9.56.dna.toplevel \
  $read_path/s_8_1_sequence.txt $read_path/s_8_2_sequence.txt

exit

# -----------------------------------------------------------------------------
# TRITUME
# -----------------------------------------------------------------------------


  --best \
  --strata \
  --seedmms 3 \

  $read_path/test_BMDM_CTRL_14JUN10_300BP_s_7_1_sequence.txt,$read_path/test_BMDM_LPS_14JUN10_300BP_s_8_1_sequence.txt $read_path/test_BMDM_CTRL_14JUN10_300BP_s_7_2_sequence.txt,$read_path/test_BMDM_LPS_14JUN10_300BP_s_8_2_sequence.txt

  -G /exports/work/vet_roslin_nextgen/dario/GTF/Sus_scrofa.Sscrofa9.56.gtf 
  $read_path/test_BMDM_CTRL_14JUN10_300BP_s_7_1_sequence.txt $read_path/test_BMDM_CTRL_14JUN10_300BP_s_7_2_sequence.txt 
