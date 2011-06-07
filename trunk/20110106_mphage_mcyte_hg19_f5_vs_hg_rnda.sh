#!/bin/bash

# qsub -P vet_roslin -pe memory 1 -l h_rt=00:10:00 job_20110106_mphage_mcyte_hg19_f5_vs_hg_rnda.sh

ref_db=/exports/work/vet_roslin_nextgen/dario/bwa/indexes/Human_ribosomal_DNA_complete_repeating_unit.fa
fastq=/exports/work/vet_roslin_nextgen/dario/bwa/output/20110106_fantom5_macro_monocyte_vs_ss9/mphage_mcyte_f5_hg19_21bp.fq
out_bwa=/exports/work/vet_roslin_nextgen/dario/bwa/output/20110106_fantom5_macro_monocyte_vs_hg_rdna/20110106_mphage_mcyte_hg19_f5_vs_hg_rnda_21bp.bwa
out_sam=/exports/work/vet_roslin_nextgen/dario/bwa/output/20110106_fantom5_macro_monocyte_vs_hg_rdna/20110106_mphage_mcyte_hg19_f5_vs_hg_rnda_21bp.sam

/exports/work/vet_roslin_nextgen/dario/bwa/bwa-0.5.8c/bwa aln $ref_db $fastq > $out_bwa

/exports/work/vet_roslin_nextgen/dario/bwa/bwa-0.5.8c/bwa samse $ref_db $out_bwa $fastq > $out_sam
	