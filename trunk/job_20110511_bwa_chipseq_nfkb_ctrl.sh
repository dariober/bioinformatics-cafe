#!/bin/bash

#
# Align chipseq reads from Mark Barnett from NFkB precipitation. See labbook 11/05/2011
#
# qsub -P vet_roslin -pe memory 2 -l h_rt=10:00:00 job_20110511_bwa_chipseq_nfkb_ctrl.sh

cd /exports/work/vet_roslin_nextgen/dario/bwa/output/20110511_mouse_chipseq_mb_nfkb

source ~/.bash_profile

# -------------------------[Input/Output]--------------------------------------

ref_genome='/exports/work/vet_roslin_nextgen/dario/bwa/indexes/Mus_musculus.mm9.ucsc.chromFa.fa'

fastq1='/exports/work/vet_roslin_nextgen/markb/Fastq/110124_EBRI093151_0065/423/con_nf_s_3_1_sequence.txt'
fastq2='/exports/work/vet_roslin_nextgen/markb/Fastq/110124_EBRI093151_0065/423/con_nf_s_3_2_sequence.txt'

fastq1_bwa='con_nf_1.bwa'
fastq2_bwa='con_nf_2.bwa'
aligned_sam='con_nf.sam'
aligned_bam='con_nf.bam'
aligned_q15_bam='con_nf_q15.bam'
aligned_q15_sorted_bam='con_nf_q15.sorted'
aligned_q15_bed='con_nf_q15.bed'

# -------------------------------[Run job]--------------------------------------

bwa aln -t 2 $ref_genome $fastq1 > $fastq1_bwa
bwa aln -t 2 $ref_genome $fastq2 > $fastq2_bwa
bwa sampe $ref_genome $fastq1_bwa $fastq2_bwa $fastq1 $fastq2 > $aligned_sam

samtools view -b -S $aligned_sam > $aligned_bam 
samtools view -q 15 -b $aligned_bam > $aligned_q15_bam
samtools sort -m 1000000000  $aligned_q15_bam $aligned_q15_sorted_bam
samtools index "${aligned_q15_sorted_bam}.bam"
# bamToBed -i "${aligned_q15_sorted_bam}.bam" > $aligned_q15_bed
## Remove *.sam and unsorted *.bam
