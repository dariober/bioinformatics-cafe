#!/usr/local/bin/python

"""
 Executes bwa and samtools on a batch of files
"""

## bsub -M 4194304 -R "rusage[mem=4096]" -n 4 -o log.txt < bwa_align.py

import os
import subprocess
import sys
import re

# ----------------------------[ Settings ]--------------------------------------

GENOME= '/lustre/reference_data/genomes/Homo_sapiens_NCBI_v36.3/hsa.fa'
FQ_DIR= '/lustre/sblab/berald01/projects/20111206_chipseq_foxm1_debbies/fastq'
FASTQ_EXT= '-strip-.*' ## Regex to get fastq files. This match will be stripped off the name and the rest used for the output files. MEMO: regex '\.s_[1-8]_sequence\.txt\.gz$' to strip full Solexa name from suffix.
BWA_ALN_OPT= ''
OUTPUT_DIR= os.path.join(os.path.split(FQ_DIR)[0], 'bwa_aligned') ## os.path.join(os.path.split(FQ_DIR)[0], 'bwa_aligned') will go to one level up the fastq dir and will create it.

SCP_DIR= 'berald01@uk-cri-lsrv04:/data01/sblab/users/berald01/bamfiles' ## At the end of the process send all the output files (.bam, .bai, .sh, .sh.log) here. Format: "username@hostname:/dest/dir". destination dir MUST exist.

genome_dict= { '/lustre/reference_data/genomes/Mus_musculus_NCBI_v37/mmu.fa': 'mm9',
               '/lustre/reference_data/genomes/Homo_sapiens_NCBI_v36.3/hsa.fa': 'hg18'} ## Tag for output files to look like '/lustre/.../fastq/bham-470.bwa.mm9.bam'

# ------------------------------------------------------------------------------

if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

fastq_files= os.listdir(FQ_DIR)
genome_tag= genome_dict[GENOME]

n= 1
for fq in fastq_files:
    """
    Prepare output file names
    """
    basename= re.sub(FASTQ_EXT, '', fq)
    fq_path= os.path.join(FQ_DIR, fq)
    sai= os.path.join(OUTPUT_DIR, fq + '.bwa.' + genome_tag + '.sai')
    sam= os.path.join(OUTPUT_DIR, re.sub(FASTQ_EXT, '.bwa.%s.sam' %genome_tag, fq) )
    bam= os.path.join(OUTPUT_DIR, re.sub(FASTQ_EXT, '.bwa.%s.unsorted.bam' %genome_tag, fq) )
    bam_sorted= os.path.join(OUTPUT_DIR, re.sub(FASTQ_EXT, '.bwa.%s' %genome_tag, fq) )
    cmd_file= os.path.join(OUTPUT_DIR, fq + '.sh')
    """
    Compile commands
    """
    bwa_aln= '/home/brown22/local/bin/bwa aln %s %s %s > %s' %(BWA_ALN_OPT, GENOME, fq_path, sai)
    bwa_samse= '/home/brown22/local/bin/bwa samse %s %s %s > %s' %(GENOME, sai, fq_path, sam)
    sam2bam= '/home/berald01/applications/samtools/samtools view -bS %s > %s' %(sam, bam)
    bam_sort= '/home/berald01/applications/samtools/samtools sort %s %s' %(bam, bam_sorted)
    bam_index= '/home/berald01/applications/samtools/samtools index %s' %(bam_sorted + '.bam')
    remove= """rm %s; rm %s; rm %s """ %(sai, sam, bam) 
    scp_cmd= 'scp %s* %s' %(os.path.join(OUTPUT_DIR, basename), SCP_DIR)
    """
    Write commands to file(s)
    """
    shfile= open(cmd_file, 'w')
    shfile.write('#!/bin/sh' + '\n\n' + 'set -e' + '\n\n')
    shfile.write(bwa_aln + '\n')
    shfile.write(bwa_samse + '\n\n')
    shfile.write(sam2bam + '\n')
    shfile.write(bam_sort + '\n')
    shfile.write(bam_index + '\n\n')
    shfile.write(remove + '\n\n')
    shfile.write(scp_cmd)
    shfile.close()
    """
    Submit jobs
    """
    bjob= 'bsub -R "rusage[mem=4096]" -J %s -n 1 -o %s < %s' %('bwa-' + str(n).zfill(2), cmd_file + '.log', cmd_file)
    n += 1
    subprocess.Popen(bjob,shell= True)