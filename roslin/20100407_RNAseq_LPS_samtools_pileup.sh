# ----------------[ 20100406_RNAseq_LPS_sscrofa9.56 ]--------------------------
# 7 April 2010
# Create a pileup file from the output of bowtie in SAM format.
# Generates also BAM and sorted BAM files
#

cd /exports/work/vet_roslin_nextgen/dario/samtools/samtools-0.1.7_x86_64-linux

## Create BAM file
./samtools view -bt /exports/work/vet_roslin_nextgen/dario/genomes/sscrofa9.56/Sscrofa9.56_Human_rDNA.fa /exports/work/vet_roslin_nextgen/dario/bowtie/output/20100406_RNAseq_LPS_sscrofa9.56/20100406_RNAseq_LPS_sscrofa9.56.sam -o /exports/work/vet_roslin_nextgen/dario/samtools/output/20100406_RNAseq_LPS_sscrofa9.56/20100406_RNAseq_LPS_sscrofa9.56.bam

## Sort BAM
./samtools sort -m 1500000000 /exports/work/vet_roslin_nextgen/dario/samtools/output/20100406_RNAseq_LPS_sscrofa9.56/20100406_RNAseq_LPS_sscrofa9.56.bam /exports/work/vet_roslin_nextgen/dario/samtools/output/20100406_RNAseq_LPS_sscrofa9.56/20100406_RNAseq_LPS_sscrofa9.56.sorted

## Pileup
./samtools pileup -s -c -f /exports/work/vet_roslin_nextgen/dario/genomes/sscrofa9.56/Sscrofa9.56_Human_rDNA.fa /exports/work/vet_roslin_nextgen/dario/samtools/output/20100406_RNAseq_LPS_sscrofa9.56/20100406_RNAseq_LPS_sscrofa9.56.sorted.bam > /exports/work/vet_roslin_nextgen/dario/samtools/output/20100406_RNAseq_LPS_sscrofa9.56/20100406_RNAseq_LPS_sscrofa9.56.pileup