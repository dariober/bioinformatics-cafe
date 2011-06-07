# -----------------------------------------------------------------------------
# Prepares a pileup file from the output of tophat (accepted_hits.sam)
# 
# Note that a SAM header is added accepted_hits.sam
#
# Tophat aligned RNAseq CTRL and LPS in separate runs (22/1/10, using GFF files)
# See Tophat run 22/1/10
#
# Tophat outputs concatenated and passed to samtools for pileup
# -----------------------------------------------------------------------------


## Path to dir samtools
samt='/exports/work/vet_roslin_nextgen/dario/samtools/samtools-0.1.7_x86_64-linux/'

## pileup SAM input (SAM alignment)
samin_ctrl='/exports/work/vet_roslin_nextgen/dario/tophat/output/20100122_RNAseq_CTRL/accepted_hits.sam'
samin_lps='/exports/work/vet_roslin_nextgen/dario/tophat/output/20100122_RNAseq_LPS/accepted_hits.sam'

## pileup output dir
samout='/exports/work/vet_roslin_nextgen/dario/samtools/output/20100426_RNAseq_CTRL-LPS_tophat_conc/'

# ------------------------------[ Start ]--------------------------------------

## Concatenate accepted_hits.sam files and add sam header (suitable for samtools)
cat /exports/work/vet_roslin_nextgen/dario/genomes/sscrofa9.56/Sscrofa9.56_Human_rDNA.samheader \
    $samin_ctrl $samin_lps > $samout/accepted_hits_header.sam

## Create BAM file
$samt/samtools view -bt /exports/work/vet_roslin_nextgen/dario/genomes/sscrofa9.56/Sscrofa9.56_Human_rDNA.fa \
   $samout/accepted_hits_header.sam \
   -o $samout/20100426_RNAseq_LPS-CTRL_conc_sscrofa9.56.bam

## Sort BAM
$samt/samtools sort -m 4000000000 \
  $samout/20100426_RNAseq_LPS-CTRL_conc_sscrofa9.56.bam \
  $samout/20100426_RNAseq_LPS-CTRL_conc_sscrofa9.56.sorted

## Pileup
$samt/samtools pileup -s -c -f \
  /exports/work/vet_roslin_nextgen/dario/genomes/sscrofa9.56/Sscrofa9.56_Human_rDNA.fa \
  $samout/20100426_RNAseq_LPS-CTRL_conc_sscrofa9.56.sorted.bam \
  > $samout/20100426_RNAseq_LPS-CTRL_conc_sscrofa9.56.pileup

## Tidy up
#rm $samout/accepted_hits_header.sam
#rm $samout/accepted_hits_header2.sam
#rm $samout/20100426_RNAseq_LPS-CTRL_conc_sscrofa9.56.bam

