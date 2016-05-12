#!/usr/bin/env python

import subprocess 
import tempfile
import os
import atexit
import sys

VERSION= '0.1.0'

docstring= """
DESCRIPTION
Intersects two bed files and return the portions overlapping:
only the a-file, only the b-file, the intersection. E.g.

|-----------|        a.bed
   |--|   |-------|  b.bed
|..|                 <- results
   |..|
       |..|
          |.|
            |-----|

Printed to stderr is the number of segments in each partition and their cumulative length

USAGE
partitionBed.py <a.bed> <b.bed>

REQUIRES
* bedtools on PATH

See also 
https://github.com/dariober/bioinformatics-cafe/tree/master/partitionBed

Version %s
""" %(VERSION)

def checkExit(process, cmd):
    stderr= p.communicate()
    if p.returncode != 0:
        sys.stderr.write(cmd + '\n')
        sys.stderr.write(str(stderr) + '\n')
        sys.exit(p.returncode)

if len(sys.argv) < 3 or sys.argv[1] == '-h':
    sys.stderr.write(docstring + '\n')
    sys.exit(1)

cmd= """subtractBed -a %s -b %s""" %(sys.argv[1], sys.argv[2])
p= subprocess.Popen(cmd, shell= True, stdout= subprocess.PIPE, stderr= subprocess.PIPE)
nseg, segLen= 0, 0
for line in iter(p.stdout.readline, b''):
    bedLine= line.strip().split('\t')
    nseg += 1
    segLen += int(bedLine[2]) - int(bedLine[1])
    sys.stdout.write(line.strip() + '\ta\n')
sys.stderr.write(str(nseg) + "\t" + str(segLen) + "\t" + sys.argv[1] + "\n")
checkExit(p, cmd)

cmd= """subtractBed -b %s -a %s""" %(sys.argv[1], sys.argv[2])
nseg, segLen= 0, 0
p= subprocess.Popen(cmd, shell= True, stdout= subprocess.PIPE, stderr= subprocess.PIPE)
for line in iter(p.stdout.readline, b''):
    bedLine= line.strip().split('\t')
    nseg += 1
    segLen += int(bedLine[2]) - int(bedLine[1])
    sys.stdout.write(line.strip() + '\tb\n')
sys.stderr.write(str(nseg) + "\t" + str(segLen) + "\t" + sys.argv[2] + "\n")
checkExit(p, cmd)

cmd= """intersectBed -a %s -b %s""" %(sys.argv[1], sys.argv[2])
nseg, segLen= 0, 0
p= subprocess.Popen(cmd, shell= True, stdout= subprocess.PIPE, stderr= subprocess.PIPE)
for line in iter(p.stdout.readline, b''):
    bedLine= line.strip().split('\t')
    nseg += 1
    segLen += int(bedLine[2]) - int(bedLine[1])
    sys.stdout.write(line.strip() + '\tab\n')
sys.stderr.write(str(nseg) + "\t" + str(segLen) + "\t" + sys.argv[1] + "," + sys.argv[2] + "\n")
checkExit(p, cmd)

