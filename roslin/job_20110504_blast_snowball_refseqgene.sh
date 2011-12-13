#!/bin/bash

# See labbook 04/05/2011
#
# Blast the sequences from which the 'Snowball' array has been design from. 
# Sequences from Chris Tuggle (non_PSR_fastas.zip, three files, concatenated into non_psr_fastas_cat.fa)
# Blasted against refseq genes from NCBI, pre-formatted blast db (ftp://ftp.ncbi.nih.gov/blast/db/refseqgene.tar.gz)
#

# qsub -P vet_roslin -l h_rt=06:00:00 job_20110504_blast_snowball_refseqgene.sh

source ~/.bash_profile
 
cd /exports/work/vet_roslin_nextgen/dario/blast/output/20110504_snowball

blast_db='/exports/work/vet_roslin_nextgen/dario/blast/db/refseqgene'
blastn -query non_psr_fastas_cat.fa -db $blast_db -task megablast -evalue 10 -out non_psr_fastas_cat.blastout -outfmt 6
## blastn -query non_psr_fastas_cat.fa -db $blast_db -task megablast -evalue 0.01 -out non_psr_fastas_cat.blastout -outfmt 6