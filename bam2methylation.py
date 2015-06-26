#!/usr/bin/env python

import sys
import argparse
import os
import subprocess
import tempfile
import shutil
import atexit

parser = argparse.ArgumentParser(description= """
DESCRIPTION
    Extract methylation calls from BAM file. Only positions with non-zero count are printed.
    If you want to get only the Cs in a certain context (e.g. CpG) you need to pass
    a bed file where these positions are. For example, to get all the
    CpGs you could use:
    fastaRegexFinder.py -f mmu.fa -r CG --noreverse > mm9.allcpg.bed
    (Possibly pipe to `cut -f1-3 | gzip` to make the file smaller)

OUTPUT:
   bedGraph with columns:
    <chrom>  <pos-1>  <pos>  <pct meth'd>  <cnt methylated>  <tot count>  <strand> [cnt mismatch]
    Memo: bedGraph is 0-based, so if the first base of chr1 is C it will have position: `chrom 0 1 ... +`
   Output is sorted by chromosome and position.
   
   With --mismatch option an additional column of mismatch count is after <tot count>.
   This is the number of A or G when the reference is C.
   The percent mismatch is given by: $7/($6+$7)

REQUIRES:
      - BAM file sorted and indexed
      - samtools on path
      - bedtools on path
      - Unix sort and awk

""", formatter_class= argparse.RawTextHelpFormatter, prog= os.path.basename(__file__))

parser.add_argument('--input', '-i',
                   required= True,
                   help='''Input bam file.
                   ''')

parser.add_argument('--ref', '-r',
                   required= True,
                   help='''Reference fasta file.
                   ''')

parser.add_argument('--l', '-l',
                   required= False,
                   help='''Bedfile with intervals where pileup should be generated.
Passed to `samtools mpileup -l <>`
                   ''')

parser.add_argument('--A', '-A',
                   action= 'store_true',
                   help='''Passed to mpileup: Count anomalous read pairs. Default
is to exclude them.
                   ''')

parser.add_argument('--samargs', '-s',
                   required= False,
                   default= '',
                   help='''String of optional arguments passed to `samtools view` to filter
reads. Put this string in quotes leaving a space after opening quote. See also --region. 
E.g. -s ' -q 15 -F 256'
                   ''')

parser.add_argument('--region', '-R',
                   required= False,
                   default= '',
                   help='''Region passed to `samtools view` to extract reads from.
E.g. -R 'chr1:0-10000'
                   ''')

parser.add_argument('--mismatch', '-mm',
                  action= 'store_true',
                  help='''Insert a column of mismatches after the tot methylation count.
                   ''')

parser.add_argument('--minq', '-mq',
                  type= int,
                  default= 0,
                  required= False,
                  help='''Minimum base quality required to consider the base a
methylation call or a mismatch.
                   ''')

parser.add_argument('--keeptmp',
                  action= 'store_true',
                  help='''Keep tmp dir. Use for debugging.
                   ''')

parser.add_argument('--version', action='version', version='%(prog)s 0.3.1')

# ------------------------------------------------------------------------------
def bam2methylation(bam, ref, bed, tmpdir, args_A, args_region, args_minq, args_mismatch, args_samargs):
   """Convert input bam to methylation files separating read 1 and read 2.
   bam:
      Input bam file, sorted and indexed
   ref:
      Reference fasta file
   bed:
      BED file of intervals to generate mpileup
   tmpdir:
      A working dir where files will be sent
   Return:
      List of length two with output file names.
   """
   if bed is None:
      L= ''
   else:
      L= '-l %s' %(bed)
   if args_A:
      A= '-A'
   else:
      A= ''
   outfilenames= []
#
# TODO: Run processes for F128 and f128 in parallel and process the outputs as they
# come through.

##   procs= []
   for F in ['-F128', '-f128']:
      ## Memo: "-F" *excludes*; "-f" *includes*
      ## -F128: Get read 1; -f128: Get read 2
      if F == '-F128':
         is_second= False
      else:
         is_second= True
      ## Prepare output
      outname= os.path.join(tmpdir, 'read%s.mpileup.txt' %(F))
      outfilenames.append(outname)
      mpileup= open(outname, 'w')
      sys.stderr.write('Methylation file: ' + outname + '\n')
      ## Prepare and execute mpileup with appropriate -F flag
      cmd_r= 'samtools view -u %(samargs)s %(F)s %(bam)s %(region)s | samtools mpileup -d100000000 -Q0 -B %(L)s -f %(ref)s %(A)s - | sort -k1,1 -s' %{
            'samargs':args_samargs, 'F': F, 'bam':bam, 'region': args_region, 'L':L, 'ref': ref, 'A': A}
      sys.stderr.write(cmd_r + '\n')
      p= subprocess.Popen(cmd_r, shell= True, stdout= subprocess.PIPE, stderr= subprocess.PIPE)
      for line in p.stdout:
          line= line.strip().split('\t')
          methList= pileup2methylation(chrom= line[0], pos= int(line[1]),
                                   callString= acceptedCalls(bases= line[4], qual_string= line[5], minq= args_minq),
                                   ref= line[2], is_second= is_second, add_mismatch= args_mismatch)

          if methList is not None:
              mpileup.write('\t'.join(methList) + '\n')
      mpileup.close()

      # Check clean exit
      stdout, stderr= p.communicate()
      if p.returncode != 0:
           print(stderr)
           print('Exit code %s' %(p.returncode))
           sys.exit(1)
   return( outfilenames )

def cleanCallString(bases):
   """Removes from the call string in mpileup (5th column) the ^ character and
   the char next to it.
   bases:
      String of read bases (5th column of mpileup)   
   Return:
      Same string as bases but with ^ and the following char removed as well as
      indels
   Example:
      bases= '^A....,,.,.,...,,,.,....^k.'
      cleanCallString(bases) >>> '....,,.,.,...,,,.,.....'
   """
   
   callString= ''
   skip= False
   getIndel= False ## Switch to start accumulating ints following +/-
   indel= []       ## List of ints following +/-. Converted to int() will give indel length
   nskip= 0
   
   for x in bases:
      if nskip > 0:
         nskip -= 1
      elif x  == '^':
         skip= True
      elif skip:
         skip= False
      elif x == '$':
         continue # Skip end-of-read marker
      elif x in ('+', '-'):
         getIndel= True
      elif getIndel:
         if x.isdigit():
            indel.append(x)
         else:
            nskip= int(''.join(indel)) - 1
            indel= []
            getIndel= False
      else:
         callString += x
   
   return(callString)
      
def pileup2methylation(chrom, pos, callString, ref, is_second= False, add_mismatch= False):
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

pileup2methylation('chr1', 1, '', C, is_second= False)
   """
   cnt_M= 0 ## Count methylated
   cnt_m= 0 ## Count unmethylated
   cnt_MM= 0 ## Count mismatches. I.e. not C or T when reference has C
   if ref.upper() == 'G':
      strand= '-'
      if is_second:
         cnt_M += callString.count('.')
         cnt_m += callString.count('A')
         cnt_MM += (callString.count('C') + callString.count('T') + callString.count('N'))
      else:
         cnt_M += callString.count(',')
         cnt_m += callString.count('a')
         cnt_MM += (callString.count('c') + callString.count('t') + callString.count('n'))
   elif ref.upper() == 'C':
      strand= '+'
      if is_second:
         cnt_M += callString.count(',')
         cnt_m += callString.count('t')
         cnt_MM += (callString.count('a') + callString.count('g') + callString.count('n'))
      else:
         cnt_M += callString.count('.')
         cnt_m += callString.count('T')
         cnt_MM += (callString.count('A') + callString.count('G') + callString.count('N'))
   else:
      return(None)
   if (cnt_m + cnt_M + cnt_MM) == 0:
      return(None)
   totreads= cnt_M + cnt_m
   
   if(totreads == 0):
      pct_met= 0
   else:
      pct_met= round(100*(float(cnt_M)/totreads), 4)
   methList= [chrom, str(pos-1), str(pos), str(pct_met), str(cnt_M), str(totreads), strand]
   
   if add_mismatch:
      methList.append(str(cnt_MM))
   
   return(methList)

def mergeMpileup(metCall_r1, metCall_r2, add_mismatch):
   """Merge the methylation call files from read 1 and read 2.
   """
   if add_mismatch:
         cmd= '''sort -m -s -k1,1 -k2,2n -k3,3n %s %s \
           | bedtools groupby -g 1,2,3 -c 5,6,7,8 -o sum,sum,distinct,sum \
           | awk '{if($4==0 && $5==0){pct= 0} else {pct= 100*($4/$5)} printf("%%s\t%%s\t%%s\t%%0.2f\t%%s\t%%s\t%%s\t%%s\\n", $1, $2, $3, pct, $4, $5, $6, $7)}'
      ''' %(metCall_r1, metCall_r2)
   else:
      cmd= '''sort -m -s -k1,1 -k2,2n -k3,3n %s %s \
           | bedtools groupby -g 1,2,3 -c 5,6,7 -o sum,sum,distinct \
           | awk '{if($4==0 && $5==0){pct= 0} else {pct= 100*($4/$5)} printf("%%s\t%%s\t%%s\t%%0.2f\t%%s\t%%s\t%%s\\n", $1, $2, $3, pct, $4, $5, $6)}'
      ''' %(metCall_r1, metCall_r2)
   sys.stderr.write(cmd + '\n')
   p= subprocess.Popen(cmd, shell= True, stdout= subprocess.PIPE, stderr= subprocess.PIPE)   
   for line in iter(p.stdout.readline, b''):
      sys.stdout.write(line)

def rmLowQualsCalls(call_string, qual_string, minq= 0):
   """Removes base calls from a string of calls if their quality is lower than
   minq.
   call_string:
      String of base calls, i.e. 5th column in mpileup. This string cleaned by
      cleanCallString() to remove mapping qualities and start/end of read markers
   qual_string:
      Quality string. i.e. 6th column from mpileup
   minq:
      Remove calls with quality lower than this. Base quality= ord(ascii) - 33
      E.g. ord(I) - 33 = 40

   Return:
      call_string with low quality bases removed.
   
   qual_string= '''!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ!"#$'''
   call_string=   'AAAAAAAAAAAAAAAAAAAATCCCCCCCCCCCCCCCCCCCCGNNNN'
   rmLowQualsCalls(call_string, qual_string, minq= 20)
   'TCCCCCCCCCCCCCCCCCCCCG'
   """
   if(len(qual_string) != len(call_string)):
      sys.exit('Quality and call string differ in length')
   
   if (minq <= 0):
      return(call_string)
   
   to_keep= []
   for i in range(0, len(qual_string)):
      q= ord(qual_string[i]) - 33
      if q >= minq:
         to_keep.append(call_string[i])
   return(''.join(to_keep))

def acceptedCalls(bases, qual_string, minq):
   """Parse the "bases" string (5th column mpileup) to return the characters to be
   used for methylation and mismatch calling
   
   bases:
      String corresponding to 5th column mpileup
   qual_string:
      Quality string (6th column mpileup) to use to filter out bad calls
   minq:
      Minimum acceptable quality to use for methylation and mismatch calling
   
   Return: String suitable for pileup2methylation() and pileup2mismatch().

bases=   '^kAAAAAAAAAAAAAAAAAAAATCCCCCCCCCCCCCCCCCCCCGNNNN^k$'   
qual_string= '''!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ!"#$'''
acceptedCalls(bases, qual_string, minq)
   'TCCCCCCCCCCCCCCCCCCCCG'
   
   """
   accepted= rmLowQualsCalls(cleanCallString(bases), qual_string, minq)
   return(accepted)
# ------------------------------------------------------------------------------

if __name__ == '__main__':

   args= parser.parse_args()
   if args.region:
      region= "'" + args.region + "'"
   else:
      region= ''
   # ------------------------------------------------------------------------------
   if not os.path.isfile(args.input):
      sys.stderr.write('\nError: File %s not found\n\n' %(args.input))
      sys.exit(1)
   if not os.path.isfile(args.ref):
      sys.stderr.write('\nError: File %s not found\n\n' %(args.ref))
      sys.exit(1)
   
   tmpdir= tempfile.mkdtemp(prefix= 'tmp_bam2methylation_')
   if not args.keeptmp:
       atexit.register(shutil.rmtree, tmpdir)
   outpiles=  bam2methylation(bam= args.input, ref= args.ref, bed= args.l, tmpdir= tmpdir, args_A= args.A, args_region= region, args_minq= args.minq, args_mismatch=args.mismatch, args_samargs= args.samargs)
   #outpiles= bam2methylation(bam= args.input, ref= args.ref, bed= args.l, tmpdir= tmpdir)
   mergeMpileup(outpiles[0], outpiles[1], args.mismatch)
   sys.exit()
