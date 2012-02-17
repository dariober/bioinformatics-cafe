#!/usr/local/bin/python

import sys
import shlex
import subprocess

if len(sys.argv) < 2:
    sys.exit("""
Given a list of bed files, calculates the pairwise overlap.
Column 4 in output answers the question: How many bases in file a are not covered by file b?
Column 5in output answers the question: What percentage of file a is covered by file b?

USAGE
    pairwise_bed_overlap.py file1.bed file2.bed file3.bed ...

    ## Get all bed files in current dir
    pairwise_bed_overlap.py `ls *.bed`

OUTPUT
    To stdout tab-separated (note that also each file with itself is returned):
file1  <bases in 1>  file1 <bases in 1>      <"bases in 1" minus "bases in 1">   <% file 1 covered by file 1>
file1  <bases in 1>  file2 <bases in 2>      <"bases in 1" minus "bases in 2">   <% file 1 covered by file 2>
file1  <bases in 1>  file3 <bases in 3>      <"bases in 1" minus "bases in 3">   <% file 1 covered by file 3>
...
fileN  <bases in N>  fileN-1 <bases in N-1>  <"bases in N" minus "bases in N-1"> <% file N covered by file N-1>

REQUIREMENTS:

sumbed.py on path (svn bioinformatics-misc)
bedtools subtractBed

It is assumed taht sumbed.py returns a string tab separated as: <filenaname> <span> <...other fileds ignored>

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
            diff= 0
            inters= filespans[fa]
        else:
            cmd= 'subtractBed -a %s -b %s | sumbed.py -' %(fa, fb)
            p= subprocess.Popen(cmd, shell= True, stdout= subprocess.PIPE)
            p.wait()
            sumout= p.stdout.read().rstrip('\n\r')
            diff= int(sumout.split('\t')[1])
            cmd= 'intersectBed -a %s -b %s | sumbed.py -' %(fa, fb)
            p= subprocess.Popen(cmd, shell= True, stdout= subprocess.PIPE)
            p.wait()
            sumout= p.stdout.read().rstrip('\n\r')
            inters= int(sumout.split('\t')[1])
            
        pair= pairwise()
        pair.filea= fa
        pair.fileb= fb
        pair.filea_span= filespans[fa]
        pair.fileb_span= filespans[fb]
        pair.diff= diff
        perc_cov= round(100*((filespans[fa] - diff) / float(filespans[fa])), 2)
        line= ([pair.filea, pair.filea_span, pair.fileb, pair.fileb_span, pair.diff, perc_cov, inters])
        line= [str(x) for x in line]
        print('\t'.join(line))
sys.exit()