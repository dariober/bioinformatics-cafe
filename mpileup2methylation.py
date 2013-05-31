#!/usr/bin/env python

import sys
import argparse

parser = argparse.ArgumentParser(description= """
DESCRIPTION
    Extract methylation calls from samtools mpileup. Only positions with non-zero count
    are printed. Output has is bedGraph with additional columns:
    <chrom>  <pos-1>  <pos>  <pct meth'd>  <cnt methylated>  <tot count>  <strand>

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

def cleanCallString(bases):
   """Removes from the call string in mpileup (5th column) the ^ character and
   the char next to it. Note TODO: Account for insertions/deletions.
   bases:
      String of read bases (5th column of mpileup)   
   Return:
      Same string as bases but with ^ and the following char removed
   Example:
      bases= '^A....,,.,.,...,,,.,....^k.'
      cleanCallString(bases) >>> '....,,.,.,...,,,.,.....'
   """
   callString= ''
   skip= False
   for x in bases:
      if x  == '^':
         skip= True
      elif skip:
         skip= False
      else:
         callString += x         
   return(callString)
      
def pileup2methylation(chrom, pos, callString, ref, outfmt):
   """Count methylated and unmethylated calls.
   chrom, pos:
      Chromosome (string) and position (int) on the pileup
   callString:
      String of bases obtained by cleanCallString
   ref:
      Reference base as obtained from 3nd column of mpileup
      
   Memo: mpileup input looks like this:
   chr7    3002089 C       2       .^~.    IA
   chr7    3002090 G       2       ..      HE
   chr7    3002114 C       2       ..      HE

   """
#   ref= mpileup[2]
#   callstring= mpileup[4]
   cnt_M= 0 ## Count methylated
   cnt_m= 0 ## Count unmethylated

   if ref.upper() == 'G':
      strand= '-'
      cnt_M += callString.count(',')
      cnt_m += callString.count('a')
   elif ref.upper() == 'C':
      strand= '+'
      cnt_M += callString.count('.')
      cnt_m += callString.count('T')
   else:
      return(None)
   if (cnt_m + cnt_M) == 0:
      return(None)
   ## Use for first version of methylation2mpileup.py
   ## methList= [mpileup[0], mpileup[1], strand, str(cnt_M), str(cnt_m), ref]
   if outfmt == 'bismark':
      methList= [chrom, str(pos), strand, str(cnt_M), str(cnt_m)]
   elif outfmt == 'bedgraph':
      totreads= cnt_M + cnt_m
      methList= [chrom, str(pos-1), str(pos), str(round(100*(float(cnt_M)/totreads), 4)), str(cnt_M), str(totreads), strand]
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
   callString= cleanCallString(line[4])
   methList= pileup2methylation(chrom= line[0], pos= int(line[1]), callString= callString, ref= line[2], outfmt= args.outfmt)
   if methList is None:
      pass
   else:
      try:
         print('\t'.join(methList))
      except IOError as e:
         break
fin.close()
sys.exit()


