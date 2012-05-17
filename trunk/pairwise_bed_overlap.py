#!/usr/local/bin/python

import sys
import shlex
import subprocess

if len(sys.argv) < 2:
    sys.exit("""
Given a list of bed files, calculates the pairwise overlap.

Columns:
    1. File A
    2. Bases spanned by A
    3. File B
    4. Bases spanned by B
    5. No. bases in A covered by B (intersection)
    6. Percentage of A covered by B (i.e. col-5 / col-2) 
    7. No. bases in A not covered by B (subtraction)
    8. Percentage of A not covered by B (i.e. col-7 / col-2)

USAGE
    pairwise_bed_overlap.py file1.bed file2.bed file3.bed ...

    ## Get all bed files in current dir
    pairwise_bed_overlap.py `ls *.bed`

REQUIREMENTS:

sumbed.py on path (svn bioinformatics-misc)
bedtools subtractBed

It is assumed that sumbed.py returns a string tab separated as: <filenaname> <span> <...other fileds ignored>

""")

class pairwise:
    def __init__(self):
        self.filea= ''
        self.filea_span= 0
        self.fileb= ''
        self.fileb_span= 0
        self.diff= 0
        
    
## files= shlex.split(sys.argv[1])
files= sys.argv[1:]
files= sorted(set(files))
## Number of bases spanned by each bed file
filespans= {}
for f in files:
    cmd= 'sumbed.py %s' %(f)
    p= subprocess.Popen(cmd, shell= True, stdout= subprocess.PIPE)
    p.wait()
    sumout= p.stdout.read().rstrip('\n\r')
    sumout= sumout.split('\t')
    filespans[sumout[0]]= int(sumout[1])

sumout= []
for fa in files:
    for fb in files:
        if fa == fb:
            continue
            #diff= 0
            #inters= filespans[fa]
        else:
            cmd= 'subtractBed -a %s -b %s | sumbed.py -' %(fa, fb)
            p= subprocess.Popen(cmd, shell= True, stdout= subprocess.PIPE)
            p.wait()
            sumout= p.stdout.read().rstrip('\n\r')
            diff= int(sumout.split('\t')[1])
#            cmd= 'intersectBed -a %s -b %s | sumbed.py -' %(fa, fb)
#            p= subprocess.Popen(cmd, shell= True, stdout= subprocess.PIPE)
#            p.wait()
#            sumout= p.stdout.read().rstrip('\n\r')
            
        pair= pairwise()
        pair.filea= fa
        pair.fileb= fb
        pair.filea_span= filespans[fa]
        pair.fileb_span= filespans[fb]
        pair.diff= diff              ## Number of bases of file a NOT covered by b
        inters= filespans[fa] - diff ## Number of bases of file a ALSO covered by b
        perc_cov= round(100*(inters / float(filespans[fa])), 2)        
        perc_unc= round(100*(diff / float(filespans[fa])), 2)
        line= ([pair.filea, pair.filea_span, pair.fileb, pair.fileb_span, inters, perc_cov, diff, perc_unc])
        line= [str(x) for x in line]
        print('\t'.join(line))
sys.exit()