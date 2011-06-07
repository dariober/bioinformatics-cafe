# -----------------------------------------------------------------------------
# Prepares a pileup file from the output of tophat (accepted_hits.sam)
# 
# Note that a SAM header is added accepted_hits.sam
#
# Tophat aligned RNAseq CTRL and LPS in the same run. See Tophat run 19/4/10
# -----------------------------------------------------------------------------

## Path to dir samtools
samt='/exports/work/vet_roslin_nextgen/dario/samtools/samtools-0.1.7_x86_64-linux/'

## pileup SAM input (SAM alignment)
samin='/exports/work/vet_roslin_nextgen/dario/tophat/output/20100419_RNAseq_LPS-CTRL/accepted_hits.sam'

## pileup output dir
samout='/exports/work/vet_roslin_nextgen/dario/samtools/output/20100426_RNAseq_CTRL-LPS_tophat/'

# ------------------------------[ Start ]--------------------------------------

## Remove SAM header from accepted_hits.sam and add custom one (suitable for samtools)
grep -v '^@' $samin > $samout/accepted_hits_header2.sam
cat /exports/work/vet_roslin_nextgen/dario/genomes/sscrofa9.56/Sscrofa9.56_Human_rDNA.samheader $samout/accepted_hits_header2.sam > $samout/accepted_hits_header.sam

## Create BAM file
$samt/samtools view -bt /exports/work/vet_roslin_nextgen/dario/genomes/sscrofa9.56/Sscrofa9.56_Human_rDNA.fa \
   $samout/accepted_hits_header.sam \
   -o $samout/20100426_RNAseq_LPS-CTRL_sscrofa9.56.bam

## Sort BAM
$samt/samtools sort -m 4000000000 \
  $samout/20100426_RNAseq_LPS-CTRL_sscrofa9.56.bam \
  $samout/20100426_RNAseq_LPS-CTRL_sscrofa9.56.sorted

## Pileup
$samt/samtools pileup -s -c -f \
  /exports/work/vet_roslin_nextgen/dario/genomes/sscrofa9.56/Sscrofa9.56_Human_rDNA.fa \
  $samout/20100426_RNAseq_LPS-CTRL_sscrofa9.56.sorted.bam \
  > $samout/20100426_RNAseq_LPS-CTRL_sscrofa9.56.pileup

## Tidy up
rm $samout/accepted_hits_header.sam
rm $samout/accepted_hits_header2.sam
rm $samout/20100426_RNAseq_LPS-CTRL_sscrofa9.56.bam

