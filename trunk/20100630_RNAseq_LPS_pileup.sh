# -----------------------------------------------------------------------------
# Prepares a pileup file from the output of tophat (accepted_hits.sam)
# 
# Note that a SAM header is added accepted_hits.sam
# -----------------------------------------------------------------------------

cd /exports/work/vet_roslin_nextgen/dario/samtools/samtools-0.1.7_x86_64-linux

## Add SAM header to accepted_hits.sam
#echo 'Removing header lines from accepted_hits_sam...'
#  sed '1,2d' /exports/work/vet_roslin_nextgen/dario/tophat/output/20100528_RNAseq_LPS/accepted_hits.sam > accepted_hits_noheader.sam
echo 'Replacing SAM header...'
  cat /exports/work/vet_roslin_nextgen/dario/genomes/sscrofa9.56/Sus_scrofa.Sscrofa9.56.dna.toplevel.samheader \
      accepted_hits_noheader.sam > accepted_hits_header.sam

## Create BAM file
./samtools view -bt /exports/work/vet_roslin_nextgen/dario/genomes/sscrofa9.56/Sus_scrofa.Sscrofa9.56.dna.toplevel.fa \
   accepted_hits_header.sam \
   -o /exports/work/vet_roslin_nextgen/dario/samtools/output/20100630_RNAseq_LPS/20100630_RNAseq_LPS.bam

## Sort BAM -m 4000000000 
./samtools sort \
  /exports/work/vet_roslin_nextgen/dario/samtools/output/20100630_RNAseq_LPS/20100630_RNAseq_LPS.bam \
  /exports/work/vet_roslin_nextgen/dario/samtools/output/20100630_RNAseq_LPS/20100630_RNAseq_LPS.sorted

## Pileup
./samtools pileup -s -c -f \
  /exports/work/vet_roslin_nextgen/dario/genomes/sscrofa9.56/Sus_scrofa.Sscrofa9.56.dna.toplevel.fa \
  /exports/work/vet_roslin_nextgen/dario/samtools/output/20100630_RNAseq_LPS/20100630_RNAseq_LPS.sorted.bam \
  > /exports/work/vet_roslin_nextgen/dario/samtools/output/20100630_RNAseq_LPS/20100630_RNAseq_LPS.pileup

# rm accepted_hits_header.sam
# rm accepted_hits_noheader.sam
# rm /exports/work/vet_roslin_nextgen/dario/samtools/output/20100630_RNAseq_LPS/20100630_RNAseq_LPS.bam