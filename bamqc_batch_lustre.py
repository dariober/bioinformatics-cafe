#!/home/berald01/.local/bin/python

"""
Execute bamqc.py on all the bam files under BAMQCDIR (searched recursivley)
Output files sent to BAMDIR.

bamqc.py is executed on lustre using bsub.
"""

import subprocess
import os

BAMQCDIR= '/lustre/sblab/berald01/bamqc' ## dir where output files will be
BAMDIR= '/lustre/sblab/berald01' ## Top dir where to serach for bams

p= subprocess.Popen('find %s' %(BAMDIR), shell= True, stdout=subprocess.PIPE)
allfiles= p.stdout.read().strip()
allfiles= allfiles.split('\n')
bamfiles= [x for x in allfiles if x.endswith('.bam')]
if not os.path.exists(BAMQCDIR):
    os.makedirs(BAMQCDIR)

for bam in bamfiles:
    bamqc_out= os.path.join(BAMQCDIR, os.path.splitext(os.path.split(bam)[1])[0] + '.bamqc.tsv')
    bamqc_log= os.path.join(BAMQCDIR, os.path.splitext(os.path.split(bam)[1])[0] + '.bamqc.log')
    cmd= 'bsub -o %s "bamqc.py --nograph --noheader -i %s -o %s"' %(bamqc_log, bam, bamqc_out) # -R "rusage[MEM=1024]"
    print(cmd)
    p= subprocess.Popen(cmd, shell= True)