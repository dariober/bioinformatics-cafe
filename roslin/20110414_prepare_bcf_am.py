#!/usr/bin/python

#
# Split BAM files into individual chromosomes and call variants
#

import os
import sys

wdir= '/exports/work/vet_roslin_nextgen/dario/samtools/output/20110414_rnaseq_mpileup'
os.chdir(wdir)

## BAM to parse
bams= ['am.bam', 'bmdm.bam']
bam_prefix= ['am', 'bmdm']

samdir= '/exports/work/vet_roslin_nextgen/dario/samtools/samtools-0.1.15/samtools'
bcftools_dir= '/exports/work/vet_roslin_nextgen/dario/samtools/samtools-0.1.15/bcftools/bcftools'

ref_genome='/exports/work/vet_roslin_nextgen/dario/genomes/sscrofa9.56/Sus_scrofa.Sscrofa9.56.dna.toplevel.fa'
chromosomes= ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', 'MT', 'X']

## Parse bam files:
snp_cmd= []
for chromosome in chromosomes:
    snp= 'snp' + chromosome + '.sh'
    fout= open(snp, 'w')
    fout.write('#!/bin/bash' + '\n')
    fout.write('cd ' + wdir + '\n')
    parsed= []
    for i in range(0, len(bams)):
        chrbam= chromosome + '.bam'
        parsedbam= 'chr' + chromosome + '.' + bams[i]
        samcmd= samdir + ' view -b ' + bams[i] + " '" + chromosome + "'" + ' > ' + parsedbam
        parsed.append(parsedbam)
        fout.write(samcmd + '\n')
    bcf= 'var.raw.' + chromosome + '.bcf'
    mpileup= samdir + ' mpileup -uf ' + ref_genome + ' ' + ' '.join(parsed) + ' | ' + bcftools_dir + ' view -bvcg - > ' + bcf
    fout.write(mpileup + '\n')
    fout.close()
    snp_cmd.append(snp)
## Submit job:
for snp in snp_cmd:
    job= 'qsub -P vet_roslin -l h_rt=01:59:00 ' + snp
    print(job)
    os.system('chmod 744 ' + snp)
    os.system(job)

## 'qsub -P vet_roslin -l h_rt=05:00:00' -pe memory 4 job_20110202_rnaseq_bmdm_ctrl.sh
# samtools mpileup -uf $ref_genome am_chr18.bam bmdm_chr18.bam | bcftools view -bvcg - > var.raw.bcf
#

#samtools index $am_ctrl_bam
#samtools index $am_lps_bam
#samtools index $bmdm_ctrl_bam
#samtools index $bmdm_lps_bam

#samtools view -b $am_ctrl_bam '18' > am_ctrl_chr18.bam
#samtools view -b $am_lps_bam '18' > am_lps_chr18.bam
#samtools view -b $bmdm_ctrl_bam '18' > bmdm_ctrl_chr18.bam
#samtools view -b $bmdm_lps_bam '18' > bmdm_lps_chr18.bam


