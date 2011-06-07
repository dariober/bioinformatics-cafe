#!/bin/bash

#
# Annotate SNP produced in /exports/work/vet_roslin_nextgen/dario/samtools/output/20110414_rnaseq_mpileup
#

cd /exports/work/vet_roslin_nextgen/dario/annovar/output/20110416_rnaseq_am_bmdm

## CTRL and LPS libraries concatenated and sent to mpileup/bcftools for snp calling:
am_cat_vcf='/exports/work/vet_roslin_nextgen/dario/samtools/output/20110414_rnaseq_mpileup/var.flt.am_cat.vcf'
bmdm_cat_vcf='/exports/work/vet_roslin_nextgen/dario/samtools/output/20110414_rnaseq_mpileup/var.flt.bmdm_cat.vcf'

## Prepare input file for annotate_variation.pl

convert2annovar.pl $am_cat_vcf -format vcf4 > var.flt.am_cat.annvar
