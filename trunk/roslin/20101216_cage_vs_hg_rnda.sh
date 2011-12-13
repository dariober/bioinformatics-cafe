#!/bin/bash

# qsub -P vet_roslin -pe memory 1 -l h_rt=00:10:00 job_20101216_cage_vs_hg_rnda.sh

ref_db=/exports/work/vet_roslin_nextgen/dario/bwa/indexes/Human_ribosomal_DNA_complete_repeating_unit.fa
## fastq=/exports/work/vet_roslin_nextgen/dario/bwa/output/20101216_cage_050810_iterate/CAGE_050810_21bp.fq
fastq=/exports/work/vet_roslin_nextgen/dario/fastq/20100805_CAGE_AM/CAGE_050810_no_adapter.fq
out_bwa=/exports/work/vet_roslin_nextgen/dario/bwa/output/20101216_CAGE_050810_hg_rdna/20101216_cage_vs_hg_rDNA_21bp.bwa
out_sam=/exports/work/vet_roslin_nextgen/dario/bwa/output/20101216_CAGE_050810_hg_rdna/20101216_cage_vs_hg_rDNA_21bp.sam

/exports/work/vet_roslin_nextgen/dario/bwa/bwa-0.5.8c/bwa aln -t 4 $ref_db $fastq > $out_bwa

/exports/work/vet_roslin_nextgen/dario/bwa/bwa-0.5.8c/bwa samse $ref_db $out_bwa $fastq > $out_sam
	