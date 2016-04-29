#!/usr/bin/env python

import subprocess 
import tempfile
import os
import atexit
import sys

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

USAGE
partitionBed.py <a.bed> <b.bed>

REQUIRES
* bedtools on PATH

See also 
https://github.com/dariober/bioinformatics-cafe/tree/master/partitionBed
"""

def checkExit(process, cmd):
    stderr= p.communicate()
    if p.returncode != 0:
        sys.stderr.write(cmd + '\n')
        sys.stderr.write(str(stderr) + '\n')
        sys.exit(p.returncode)


if sys.argv[1] == '-h' or len(sys.argv) != 3:
    sys.stderr.write(docstring + '\n')
    sys.exit(1)

cmd= """subtractBed -a %s -b %s""" %(sys.argv[1], sys.argv[2])
p= subprocess.Popen(cmd, shell= True, stdout= subprocess.PIPE, stderr= subprocess.PIPE)
for line in iter(p.stdout.readline, b''):
    sys.stdout.write(line.strip() + '\ta\n')
checkExit(p, cmd)

cmd= """subtractBed -b %s -a %s""" %(sys.argv[1], sys.argv[2])
p= subprocess.Popen(cmd, shell= True, stdout= subprocess.PIPE, stderr= subprocess.PIPE)
for line in iter(p.stdout.readline, b''):
    sys.stdout.write(line.strip() + '\tb\n')
checkExit(p, cmd)

cmd= """intersectBed -a %s -b %s""" %(sys.argv[1], sys.argv[2])
p= subprocess.Popen(cmd, shell= True, stdout= subprocess.PIPE, stderr= subprocess.PIPE)
for line in iter(p.stdout.readline, b''):
    sys.stdout.write(line.strip() + '\tab\n')
checkExit(p, cmd)

