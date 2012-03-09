#!/home/berald01/.local/bin/python

"""
Snippet to launch a number of shell subprocesses using a limited number
if cores.
"""

import subprocess
import os
import time

NCORES= 6 ## Max number of cores to use 

bams= [x for x in sorted(os.listdir(os.getcwd())) if x.endswith('.bam')]

procs= []
running_procs= 0
n= 0
for bam in bams:
    n += 1
    cmd= 'bamqc.py -i %s --noheader --nograph >> bamqc.tsv' %(bam)
    print('%s. %s' %(n, cmd))
    p= subprocess.Popen(cmd, shell= True)
    procs.append(p)
    running_procs += 1
    while running_procs >= NCORES:
        time.sleep(5)
        running_procs= 0
        for x in procs:
            xv= x.poll()
            if x.returncode is None:
                running_procs += 1