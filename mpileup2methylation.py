#!/usr/bin/env python

import sys
import argparse

parser = argparse.ArgumentParser(description= """
DESCRIPTION
    Extract methylation calls from samtools mpileup. Only positions with non-zero count
    are printed. Output has is bedGraph with additional columns:
    <chrom>  <pos-1>  <pos>  <pct meth'd>  <tot count>  <strand>

   Memo: bedGraph is 0-based, so if the first base of chr1 is C it will have position: `chrom 0 1 ... +`

EXAMPLE & RECIPES
    samtools mpileup -d100000 -Q0 -B -l mm9.allcpg.bed.gz -f mmu.fa myreads.bam | mpileup2methylation.py -i -
   
    If you want to get only the Cs in a certain context (e.g. CpG) you need to pass
    a bed file to mpileup where these positions are. For example, to get all the
    CpGs you could use:
    fastaRegexFinder.py -f mmu.fa -r CG --noreverse > mm9.allcpg.bed
    (Possibly pipe to `cut -f1-3 | gzip` to make the file smaller)
   
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

parser.add_argument('--outfmt', '-f',
                   required= False,
                   choices= ['bedgraph', 'bismark'],
                   default= 'bedgraph',
                   help='''Output format: bedgraph (default) or bismark.
'bismark' output has columns: chrom, pos, strand, count meth'd, count unmeth'd.
The first 5 column of the genome-wide cytosine report produced by
bismark_methylation_extractor should be nearly identical to 'bismark' format.
                   ''')

args= parser.parse_args()
# ------------------------------------------------------------------------------

def pileup2methylation(mpileup, outfmt):
   """Convert a row from mpileup (as list) to a methylation call.
   mpileup= ['chr7', '3002089', 'C', '2', '.^~.', 'IA']
   pileup2methylation(mpileup)
   
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
   ## Use for first version of methylation2mpileup.py
   ## methList= [mpileup[0], mpileup[1], strand, str(cnt_M), str(cnt_m), ref]
   if outfmt == 'bismark':
      methList= [mpileup[0], mpileup[1], strand, str(cnt_M), str(cnt_m)]
   elif outfmt == 'bedgraph':
      totreads= cnt_M + cnt_m
      methList= [mpileup[0], str(int(mpileup[1])-1), mpileup[1], str(round(100*(float(cnt_M)/totreads), 4)), str(totreads), strand]
   else:
      sys.exit('Unexpected keyword for outfmt "%s"' %(outfmt))
   return(methList)
# ------------------------------------------------------------------------------

if args.input == '-':
   fin= sys.stdin
else:
   fin= open(args.input)

for line in fin:
   line= line.strip().split('\t')
   methList= pileup2methylation(line, args.outfmt)
   if methList is None:
      pass
   else:
      try:
         print('\t'.join(methList))
      except IOError as e:
         break
fin.close()
sys.exit()