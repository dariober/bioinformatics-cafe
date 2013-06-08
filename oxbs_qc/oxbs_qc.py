#!/usr/bin/env python

import argparse
import sys
import os
import subprocess
import tempfile
import re
import gzip
import shlex
import shutil
import oxbs_qc_func

parser = argparse.ArgumentParser(description= """
DESCRIPTION
    Execute QC of BS/oxBS libraries.

PIPELINE
    - Shorten reads oxbs_qc.trim_fastq()
    - Remove adapters trim_galore
    - Map reads w/ bismark
    - Convert sam to bam, sort & index
    - Clena names & Clip overallaping reads `bam clipOverlap`
    - Call methylation 
    - Plot coverage
    - 

REQUIRES
    - fastx_trimmer
    - cutadapt
    - trim_galore
    - Bismark
    - Bowtie2
    - samtools
    - cleanBamReadNames.py -b -i
    - bamUtil clipOverlap
    - mpileup2methylation.py
    
bismark --bowtie2 -o bismark_out /lustre/sblab/berald01/reference_data/genomes/bsseq_synthetic4/ fastq_trimmed/$trimmedfq &&
mv bismark_out/${trimmedfq}_bt2_bismark.sam bismark_out/${bname}.bs4.sam &&
sam2bam.py bismark_out/${bname}.bs4.sam &&
rm fastq_trimmed/$trimmedfq

    """, formatter_class= argparse.RawTextHelpFormatter, version= '0.1.0')

# -----------------------------------------------------------------------------
input_args= parser.add_argument_group('Input options', '')

input_args.add_argument('--fastq', '-f',
                   required= True,
                   nargs= '+',
                   help='''One or two fastq files to QC. If two files are passed
they will be treated as paired-ends
                   ''')

input_args.add_argument('--refdir', '-r',
                   required= True,
                   help='''The path to the directory containing the reference
sequences. This directory is expected to contain
- The reference FASTA file <ref>.fa
- The tab delimited file of base modifications <ref>.txt
- The subdir `Bisulfite_Genome` prepared by `bismark_genome_preparation --bowtie2`
                   ''')

input_args.add_argument('--outdir', '-o',
                   required= False,
                   default= 'oxbs_qc_out',
                   help='''Directory where output will be sent. It will be created
if it doesn't exist.
                   ''')

input_args.add_argument('--skip_shorten', '-Ss',
                   action= 'store_true',
                   help='''Do not perform shortening of reads. Use this option
if your reads are 3 or more bases shorter than the shortest reference sequence.
                   ''')

input_args.add_argument('--skip_trim', '-St',
                   action= 'store_true',
                   help='''Do not perform trim reads by removing 3'-adapters and
low quality ends.
                   ''')

input_args.add_argument('--bin',
                   required= False,
                   default= '',
                   help='''Path to the executable files required by oxbs_qc. Default 
is '' (all executables are on $PATH)
                   ''')

def main():
    args = parser.parse_args()
    wdir= get_wdir(args.outdir)
    print()

#    encoding= get_fastq_encoding(args.fastq[0])
#    refsize= getReferenceSize(args.refdir)
#    fqkit= FastqKit()
#    fqkit.arg_fq= args.fastq
#    fqkit.dissect_fastqname()
#    if not args.skip_shorten:
#        print('Shortening reads: %s' %(', '.join(fqkit.arg_fq)))
#        fqkit= task_trim_fastq(fqkit, refsize, wdir)
#        print('Shortened reads in: %s' %(', '.join(fqkit.arg_fq)))
#    ## trim galore
#    if not args.skip_trim:
#        fqkit= task_trim_galore(fqkit, opts= '-o %(wdir)s' %{'wdir': wdir}, path= args.bin)
#    ## Bismark alignment:
    
        
if __name__ == '__main__':
    main()
    sys.exit()
    