## Move to a convenient dir.
cd /exports/work/vet_roslin_nextgen/dario/Tritume

## Download dataset
wget ftp://ftp.ncbi.nih.gov/pub/geo/DATA/supplementary/samples/GSM611nnn/GSM611116/GSM611116%5F2060%5Fs%5F4%5Fexport%2Etxt%2Egz

## Uncompress sequences
gunzip GSM611116_2060_s_4_export.txt.gz

#
#  Execute to reconstruct the FASTQ file.
./gsm_to_fastq.py

# -----------------------------------
# Download and index mouse genome 
# -----------------------------------

## Download mouse geneome from ENSEMBL
cd /exports/work/vet_roslin_nextgen/dario/bwa/indexes 
wget ftp://ftp.ensembl.org/pub/current/fasta/mus_musculus/dna/Mus_musculus.NCBIM37.60.dna.toplevel.fa.gz

## Uncompress
gunzip /exports/work/vet_roslin_nextgen/dario/bwa/indexes/Mus_musculus.NCBIM37.60.dna.toplevel.fa.gz

## Submit job to eddie to make index
qlogin -P vet_roslin

## Make BWA indexes
cd /exports/work/vet_roslin_nextgen/dario/bwa/indexes
/exports/work/vet_roslin_nextgen/dario/bwa/bwa-0.5.8c/bwa index -a bwtsw /exports/work/vet_roslin_nextgen/dario/bwa/indexes/Mus_musculus.NCBIM37.60.dna.toplevel.fa

# Exit qlogin when done
# exit

# -------------------------------------
# Align sequences
# -------------------------------------

# see 
# /exports/work/vet_roslin_nextgen/dario/bwa/output/20110121_mouse_chipseq_mb/job_20110121_test_aln.sh

## Have a look to alignment stats
samtools flagstat GSM611116_2060_s_4_export.bam > flagstats.txt

samtools view -q 15 -F 4 GSM611116_2060_s_4_export.bam | wc


# --------------------------------------
# Convert SAM output to BED
# --------------------------------------

## Bedtools and SAMtools
PATH=$PATH:/exports/work/vet_roslin_nextgen/dario/BEDTools/BEDTools-Version-2.10.1/bin
PATH=$PATH:/exports/work/vet_roslin_nextgen/dario/samtools/samtools-0.1.10

## If necessary prepare the .fai header with
# samtools faidx /exports/work/vet_roslin_nextgen/dario/bwa/indexes/Mus_musculus.NCBIM37.60.dna.toplevel.fa

cd /exports/work/vet_roslin_nextgen/dario/bwa/output/20110121_mouse_chipseq_mb

## Convert SAM to BAM using the appropriate .fai header:

samtools view -b -S \
    -o GSM611116_2060_s_4_export.bam \
    -t /exports/work/vet_roslin_nextgen/dario/bwa/indexes/Mus_musculus.NCBIM37.60.dna.toplevel.fa.fai \
    GSM611116_2060_s_4_export.sam

## Convert BAM to BED:

bamToBed -i GSM611116_2060_s_4_export.bam > GSM611116_2060_s_4_export.bed

# ---------------------------------------
# HOMER
# ---------------------------------------

# PATH to homer: 
# /exports/work/vet_roslin_nextgen/dario/homer/bin/

cd /exports/work/vet_roslin_nextgen/dario/homer/
mkdir output
cd output
mkdir 20110121_mouse_chipseq_mb
cd 20110121_mouse_chipseq_mb


# Make tag directory
makeTagDirectory tag_dir  \
    /exports/work/vet_roslin_nextgen/dario/bwa/output/20110121_mouse_chipseq_mb/GSM611116_2060_s_4_export.bed \
    -format bed


# Find peaks
cd /exports/work/vet_roslin_nextgen/dario/homer/output/20110121_mouse_chipseq_mb

findPeaks tag_dir -center -o GSM611116_2060_s_4_export.peaks

# Annotate peaks
## Add mm8 to available genomes
# perl /exports/work/vet_roslin_nextgen/dario/homer/.//configureHomer.pl -install /exports/work/vet_roslin_nextgen/dario/homer/data/genomes/mm8

cd /exports/work/vet_roslin_nextgen/dario/homer/output/20110121_mouse_chipseq_mb

annotatePeaks.pl GSM611116_2060_s_4_export.peaks mm8 > GSM611116_2060_s_4_export.ann