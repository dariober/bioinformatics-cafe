#!/usr/bin/python

import os

base= '/exports/work/vet_roslin_nextgen/dario/tophat/output'
dirs= ['20110202_rnaseq_am', '20110202_rnaseq_am_ctrl', '20110202_rnaseq_am_lps', '20110202_rnaseq_bmdm', '20110202_rnaseq_bmdm_ctrl', '20110202_rnaseq_bmdm_lps']

samtools= '/exports/work/vet_roslin_nextgen/dario/samtools/samtools-0.1.10/samtools '

for d in dirs:
    print('Moving to ' + d)
    cdir= os.path.join(base, d)
    os.chdir(cdir)
    fout= open('samtools_flagstat.sh', 'w')
    fout.write('#!/bin/bash\n\n')
    fout.write('cd ' + cdir + '\n\n')
    fout.write(samtools + ' flagstat accepted_hits.bam > samtools_flagstat.txt\n')
    fout.close()
    os.system('chmod 744 samtools_flagstat.sh')
    os.system('qsub -P vet_roslin samtools_flagstat.sh')