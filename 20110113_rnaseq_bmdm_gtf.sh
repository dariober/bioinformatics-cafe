#!/bin/bash

#
# Align RNAseq libraries rnaseq_20100614 lps and crl against Sscrofa genome 9 
# Unfiltered reads (OK for s_7_1, s_7_2, s_8_1; Not for s_8_2).
# Use GTF file 9.60 as refernce.
#
# Note the use of --best --strata --seedmms (not allowed in original script tophat)
#
# This command sent to eddie:
# cd /exports/work/vet_roslin_nextgen/dario/tophat/output/20101116_RNAseq_combined
# qsub -P vet_roslin -l h_rt=23:00:00 -pe memory 4 job_20110113_rnaseq_bmdm_gtf.sh

. /etc/profile
module add python/2.6.3

#----------------------[ Add bowtie and tophat to PATH ]-----------------------
# Uncomment as necessary

PATH=$PATH:/exports/work/vet_roslin_nextgen/dario/tophat/tophat-1.1.4.Linux_x86_64
PATH=$PATH:/exports/work/vet_roslin_nextgen/dario/samtools/samtools-0.1.10
PATH=$PATH:/exports/work/vet_roslin_nextgen/dario/bowtie/bowtie-0.12.7

## Where reads are
cd /exports/work/vet_roslin_nextgen/dario/fastq/20100614_RNAseq_pig_BMDM/478/


# --mate-inner-dist was calculated as [300 - (32 + 33 + 57 + 57)]

tophat \
  --output-dir /exports/work/vet_roslin_nextgen/dario/tophat/output/20110113_rnaseq_bmdm_gtf \
  --mate-inner-dist 121 \
  --mate-std-dev 50 \
  --solexa1.3-quals \
  --num-threads 4 \
  --best \
  --strata \
  --seedmms 3 \
  --library-type fr-unstranded\
  -G /exports/work/vet_roslin_nextgen/dario/ensembl/release-60/gtf/sus_scrofa/Sus_scrofa.Sscrofa9.60.gtf \
  /exports/work/vet_roslin_nextgen/dario/bowtie/indexes/Sus_scrofa.Sscrofa9.56.dna.toplevel \
  BMDM_CTRL_14JUN10_300BP_s_7_1_sequence.txt,BMDM_LPS_14JUN10_300BP_s_8_1_sequence.txt BMDM_CTRL_14JUN10_300BP_s_7_2_sequence.txt,BMDM_LPS_14JUN10_300BP_s_8_2_sequence.txt


#---------------------------[ Convert BAM to SAM with samtools ]----------------

PATH=$PATH:/exports/work/vet_roslin_nextgen/dario/samtools/samtools-0.1.10
cd /exports/work/vet_roslin_nextgen/dario/tophat/output/20110113_rnaseq_bmdm_gtf

# Parse accepted_hits.bam to extract CTRL and LPS reads.

## LPS
samtools view -h accepted_hits.bam | grep -E '^@|^EBRI093151_0059:8:' > accepted_hits_lps.sam
samtools view -S -b -h -o accepted_hits_lps.bam accepted_hits_lps.sam

## CTRL
samtools view -h accepted_hits.bam | grep -E '^@|^EBRI093151_0059:7:' > accepted_hits_ctrl.sam
samtools view -S -b -h -o accepted_hits_ctrl.bam accepted_hits_ctrl.sam

exit

#-----------------------------[ Sort by qname, rname, pos ]---------------------

. /etc/profile
module add python/2.6.3
cd /exports/work/vet_roslin_nextgen/dario/tophat/output/20101116_RNAseq_combined

python /exports/work/vet_roslin_nextgen/dario/python_scripts/sort_file-2.6.py -b 4000000 -k "(line.split('\t')[0], line.split('\t')[2], int(line.split('\t')[3]))" accepted_hits.sam accepted_hits.sam.sorted

#---------------------------[ Collect unmapped reads ]--------------------------

. /etc/profile
module add python/2.6.3
cd /exports/work/vet_roslin_nextgen/dario/tophat/output/20101116_RNAseq_combined

python fastq_unmapped.py










