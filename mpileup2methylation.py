#!/usr/bin/env python

import pysam
import sys
import argparse
import re

"""samtools mpileup -l ../../Tritume/mm9.allcpg.bed -Bf /lustre/sblab/berald01/reference_data/genomes/Mus_musculus_NCBI_v37/mmu.fa - > ../../Tritume/mjb041_E14BSAD06.mm9.cpg.pileup """


parser = argparse.ArgumentParser(description= """
DESCRIPTION
    Extract methylation calls from samtools mpileup. Only positions with non-zero count
    (either methylated or not) are printed.
    
    If you want to get only the Cs in a certain context (e.g. CpG) you need to pass
    a bed file to mpileup where these positions can be found. For example,
    to get all the CpGs you could use:
    fastaRegexFinder.py -f mmu.fa -r CG --noreverse > mm9.allcpg.bed
    (Possibly pipe to cut -f1-3 to make the file smaller)
    
EXAMPLE
    samtools mpileup -d100000 -Q0 -B -l mm9.allcpg.bed -f mmu.fa myreads.bam | mpileup2methylation.py -i -

TESTED:
   14/02/2013: The output of mpileup2methylation.py is almost (but not) identical
   to bismark_methylation_extractor (files *.genome-wide_CX|CpG_report.txt). Tested
   on synthetic strands, single end, and paired-end. Differences not investigated
   further but should affect ~1 call in 1000.

""", formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('--input', '-i',
                   required= True,
                   help='''Input file as produced by samtools mpileup. Use - to read stream from stdin.
                   ''')

args= parser.parse_args()
# ------------------------------------------------------------------------------

def pileup2methylation(mpileup):
   """Convert a row from mpileup (as list) to a methylation call.
   If the reference base is C or G, return a list as :
      <chr> <pos> <strand +|-> <count methy'd> <count unmethy'd> <ref base C|G>
   If reference base is not C or G return None
   
   Memo: mpileup input looks like this:
   chr7    3002089 C       2       .^~.    IA
   chr7    3002090 G       2       ..      HE
   chr7    3002114 C       2       ..      HE

   """
   ref= mpileup[2]
   callstring= mpileup[4]
   cnt_M= 0
   cnt_m= 0
   if ref.upper() == 'G':
      strand= '-'
      cnt_M += callstring.count(',')
      cnt_m += callstring.count('a')
   elif ref.upper() == 'C':
      strand= '+'
      cnt_M += callstring.count('.')
      cnt_m += callstring.count('T')
   else:
      return(None)
   if (cnt_m + cnt_M) == 0:
      return(None)
   methList= [mpileup[0], mpileup[1], strand, str(cnt_M), str(cnt_m), ref]
   return(methList)
# ------------------------------------------------------------------------------

if args.input == '-':
   fin= sys.stdin
else:
   fin= open(args.input)

for line in fin:
   line= line.strip().split('\t')
   methList= pileup2methylation(line)
   if methList is None:
      pass
   else:
      try:
         print('\t'.join(methList))
      except IOError as e:
         break
fin.close()
sys.exit()