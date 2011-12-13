#
#  Align RNAseq_CTRL and RNAseq_LPS libs 
#  s_8_1_sequence.txt is LPS  
#

cd /exports/work/vet_roslin_nextgen/dario/bowtie/output/

/exports/work/vet_roslin_nextgen/dario/bowtie/current/bowtie \
  /exports/work/vet_roslin_nextgen/dario/bowtie/indexes/pig_repbase_20090604_renamed \
  -1 /exports/work/vet_roslin_nextgen/dario/fastq/20091217_RNAseq_AM_LPS_CTRL/s_8_1_sequence.txt \
  -2 /exports/work/vet_roslin_nextgen/dario/fastq/20091217_RNAseq_AM_LPS_CTRL/s_8_2_sequence.txt \
  -q \
  --solexa1.3-quals \
  -a \
  --seedmms 3 \
  --best \
  --strata \
  --sam \
  /exports/work/vet_roslin_nextgen/dario/bowtie/output/20100629_RNAseq_CTRL_repbase/20100629_RNAseq_CTRL_repbase.sam