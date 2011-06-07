#!/usr/bin/python

#
#  Submit BWA jobs to eddie to align Blackface DGE to Btau4.0
#

# ---------------------------------------------------

## Where output files will be:
cwd= '/exports/work/vet_roslin_nextgen/dario/bwa/output/20101222_blackface_dge'

## Reference sequence, indexed for BWA:
ref_db= '/exports/work/vet_roslin_nextgen/dario/bwa/indexes/Bos_taurus.Btau_4.0.60.dna.toplevel.fa'

## Where FASTQ files are, spliced and converted to Sanger quality. Extension has to be '.fq'
fastq_dir= '/exports/work/vet_roslin_nextgen/dario/fastq/20101223_Blackface'

# ---------------------------------------------------

import os
os.chdir(cwd)

fastq= os.listdir(fastq_dir)
fastq= [fq for fq in fastq if fq.endswith('.fq')]

for file in fastq:
    job_name= file + '.sh'
    fastq_file= os.path.join(fastq_dir, file)
    out_bwa= file + '.bwa'
    out_sam= file + '.sam'
    job_file= open(job_name, 'w')
    job_file.write('#!/bin/bash\n\n')
    job_file.write('PATH=$PATH:/exports/work/vet_roslin_nextgen/dario/bwa/bwa-0.5.8c\n\n')
    job_file.write('cd ' + cwd + '\n\n')
    job_file.write('bwa aln -t 2 ' + ref_db + ' ' + fastq_file + ' > ' + out_bwa + '\n')
    job_file.write('bwa samse ' + ref_db + ' ' + out_bwa + ' ' + fastq_file + ' > ' + out_sam + '\n')
    job_file.close()
    os.system('chmod 744 ' + job_name)
    qsub_line= 'qsub -P vet_roslin -pe memory 2 -l h_rt=05:00:00 ' + job_name
    os.system(qsub_line)
