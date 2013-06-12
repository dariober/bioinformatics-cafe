#!/usr/bin/env python

import sys
import subprocess
import os
import argparse
import string

parser = argparse.ArgumentParser(description= """
DESCRIPTION
    Convert a SAM file to BAM, sort and index.
    NB: Input *.sam is replaced by *.bam
    
    This is what is executed:
    -------------------------
    samtools view -S -u %(sam)s | samtools sort - %(bname)s &&
    samtools index %(bname)s.bam &&
    rm %(sam)s
    
USAGE
    sam2bam.py <input.sam>
    
REQUIRES
    samtools on PATH

""", formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('input',
                   help='''Input file as produced by samtools mpileup. Use - to read from stdin.

                   ''')

parser.add_argument('--noidx',
                   action= 'store_true',
                   help='''Do not index the sorted output bam.

                   ''')

parser.add_argument('--force', '-f',
                   action= 'store_true',
                   help='''Overwrite bam files if it exists. Default is to fail in such
case.

                   ''')

parser.add_argument('--markdup', '-m',
                   action= 'store_true',
                   help='''Mark duplicates. Requires picard MarkDuplicates.jar on the PATH.
The metrics file produced by picard is written to <input w/o .sam>.markDuplicates.txt

                   ''')

parser.add_argument('--echo', '-e',
                   action= 'store_true',
                   help='''Only print the command that would be executed and exit.

                   ''')

# -----------------------------------------------------------------------------

def search_file(filename, search_path):
   """Given a search path, find file.
   See http://code.activestate.com/recipes/52224-find-a-file-given-a-search-path/
   NB: This fun could be replaced by distutils.spawn.find_executable(cmd) 
   """
   file_found = 0
   paths = string.split(search_path, os.pathsep)
   for path in paths:
      if os.path.exists(os.path.join(path, filename)):
          file_found = 1
          break
   if file_found:
      return(os.path.abspath(os.path.join(path, filename)))
   else:
      return(None)

def main():
   args= parser.parse_args()
   
   sam= args.input
   if sam == '-':
       sam= sys.stdin.readlines()
       if len(sam) > 1:
           sys.exit('Only one input file allowed!. Got %s' %(sam))
       sam= sam[0].strip()
   if not sam.endswith('.sam'):
       sys.exit('Input file "%s" does not have extension .sam' %(sam))
   
   ## Create output file names:
   bname= os.path.splitext(sam)[0]
   unsortedbam= bname + ''
   bam= bname + '.bam'
   
   if os.path.exists(bam) and not args.force:
       sys.exit('Bam file "%s" already exists. Use -f to overwrite' %(bam))
   
   ## NB: could use -u option but it throws and error (bug)
   sortedBam= bname + '.bam'
   cmd_bam= """\nsamtools view -S -b %(sam)s | samtools sort - %(bname)s &&""" %{'sam': sam, 'bname': bname}
   cmd_idx= """samtools index %(sortedBam)s &&"""  %{'sortedBam': sortedBam}
   cmd_rm= """rm %(sam)s\n""" %{'sam': sam}
   
   # Picard 
   # -----------------------------------------------------------------------------
   mdBam= bname + '.tmp.bam'
   metricsFile= bname + '.markDuplicates.txt'
   
   if args.markdup:
       markdup= search_file('MarkDuplicates.jar', os.environ['PATH'])
       if markdup is None:
           sys.exit('Cannot find MarkDuplicates.jar on your PATH')
       cmd_md= '\njava -Xmx2g -jar ' + markdup + ' INPUT="' + sortedBam + '" OUTPUT="' + mdBam + '" METRICS_FILE="' + metricsFile + '" ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT && \n' + 'mv "' + mdBam + '" "' + sortedBam + '" &&'
   else:
       cmd_md= ''
   cmd_bam += cmd_md
   # -----------------------------------------------------------------------------
   
   if args.noidx:
       cmd= '\n'.join([cmd_bam, cmd_rm])
   else:
       cmd= '\n'.join([cmd_bam, cmd_idx, cmd_rm])
   
   print(cmd)
   if not args.echo:
       p= subprocess.Popen(cmd, stdout= subprocess.PIPE, stderr= subprocess.PIPE, shell= True)
       p.wait()
       print(p.stderr.read())
   sys.exit()

if __name__ == '__main__':
   main()
   
