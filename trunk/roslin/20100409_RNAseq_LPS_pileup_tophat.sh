# -----------------------------------------------------------------------------
# Prepares a pileup file from the output of tophat (accepted_hits.sam)
# 
# Note that a SAM header is added accepted_hits.sam
# -----------------------------------------------------------------------------

cd /exports/work/vet_roslin_nextgen/dario/samtools/samtools-0.1.7_x86_64-linux

## Add SAM header to accepted_hits.sam
cat /exports/work/vet_roslin_nextgen/dario/genomes/sscrofa9.56/Sscrofa9.56_Human_rDNA.samheader \
  /exports/work/vet_roslin_nextgen/dario/tophat/output/20100122_RNAseq_LPS/accepted_hits.sam > tritume/accepted_hits_header.sam

## Create BAM file
./samtools view -bt /exports/work/vet_roslin_nextgen/dario/genomes/sscrofa9.56/Sscrofa9.56_Human_rDNA.fa \
   tritume/accepted_hits_header.sam \
   -o /exports/work/vet_roslin_nextgen/dario/samtools/output/20100409_RNAseq_LPS_tophat/20100409_RNAseq_LPS_sscrofa9.56.bam

## Sort BAM
./samtools sort -m 4000000000 \
  /exports/work/vet_roslin_nextgen/dario/samtools/output/20100409_RNAseq_LPS_tophat/20100409_RNAseq_LPS_sscrofa9.56.bam \
  /exports/work/vet_roslin_nextgen/dario/samtools/output/20100409_RNAseq_LPS_tophat/20100409_RNAseq_LPS_sscrofa9.56.sorted

## Pileup
./samtools pileup -s -c -f \
  /exports/work/vet_roslin_nextgen/dario/genomes/sscrofa9.56/Sscrofa9.56_Human_rDNA.fa \
  /exports/work/vet_roslin_nextgen/dario/samtools/output/20100409_RNAseq_LPS_tophat/20100409_RNAseq_LPS_sscrofa9.56.sorted.bam \
  > /exports/work/vet_roslin_nextgen/dario/samtools/output/20100409_RNAseq_LPS_tophat/20100409_RNAseq_LPS_sscrofa9.56.pileup

rm tritume/accepted_hits_header.sam