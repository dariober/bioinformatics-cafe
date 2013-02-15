#!/usr/bin/env python

import sys
import argparse

parser = argparse.ArgumentParser(description= """
DESCRIPTION
    Extract methylation calls from samtools mpileup. Only positions with non-zero count
    are printed. Output has 1-based coordinates with columns:
    <chrom> <pos> <strand> <count meth> <count unmeth> <C|G>
    
EXAMPLE & RECIPES
    samtools mpileup -d100000 -Q0 -B -l mm9.allcpg.bed.gz -f mmu.fa myreads.bam | mpileup2methylation.py -i -
    chr7	3002089	+	2	0	C
    chr7	3002114	+	2	0	C
    chr7	3002122	+	2	0	C
   
    If you want to get only the Cs in a certain context (e.g. CpG) you need to pass
    a bed file to mpileup where these positions are. For example, to get all the
    CpGs you could use:
    fastaRegexFinder.py -f mmu.fa -r CG --noreverse > mm9.allcpg.bed
    (Possibly pipe to `cut -f1-3 | gzip` to make the file smaller)

   Convert output of mpileup2methylation.py to bed:
   
   ... | mpileup2methylation.py -i - | awk '{print $1 "\t" $2-1 "\t" $2 "\t" $6 "\t.\t" $3 "\t" $4 "\t" $5}'
   
   chr7    3002088 3002089 C       .       +       2       0
   chr7    3002113 3002114 C       .       +       2       0
   chr7    3002121 3002122 C       .       +       2       0

   Note that the oiginal position is decreased by 1. This bed passed to fastaFromBed will return the same base as
   in the pileup.
   
TESTED:
   14/02/2013: The output of mpileup2methylation.py is almost (but not) identical
   to bismark_methylation_extractor (files *.genome-wide_CX|CpG_report.txt). Tested
   on synthetic strands single end, and paired-end whole genome. Differences not investigated
   further but should affect ~1 call in 1000-10000.

""", formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('--input', '-i',
                   required= True,
                   help='''Input file as produced by samtools mpileup. Use - to read from stdin.
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
   cnt_M= 0 ## Count methylated
   cnt_m= 0 ## Count unmethylated
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