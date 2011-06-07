#!/usr/bin/python

#
# Prepares cufflinks scripts and submits to eddie
#

import os

target_dirs= ['/exports/work/vet_roslin_nextgen/dario/tophat/output/20110202_rnaseq_am_lps/cufflinks_denovo',
              '/exports/work/vet_roslin_nextgen/dario/tophat/output/20110202_rnaseq_bmdm_ctrl/cufflinks_denovo',
              '/exports/work/vet_roslin_nextgen/dario/tophat/output/20110202_rnaseq_bmdm_lps/cufflinks_denovo']

for d in target_dirs:
    try:
        os.mkdir(d)
    except:
        pass
    cuff_fname= os.path.join(d, 'job_20110428_cufflinks_denovo.sh')
    fout= open(cuff_fname, 'w')
    
    cuff_script="""
#!/bin/bash

# -----------------------------------------------------------------------------
# Execute cufflinks for no denovo assembly
# -----------------------------------------------------------------------------

# qsub -P vet_roslin -l h_rt=02:00:00 job_20110428_cufflinks_denovo.sh

PATH=$PATH:/exports/work/vet_roslin_nextgen/dario/cufflinks/cufflinks-0.9.3.Linux_x86_64

cd %s

cufflinks ../accepted_hits.bam -r /exports/work/vet_roslin_nextgen/dario/bowtie/indexes/Sus_scrofa.Sscrofa9.56.dna.toplevel.fa -N

mv transcripts.expr transcripts_denovo.expr
mv transcripts.gtf transcripts_denovo.gtf
""" %(d)
    fout.write(cuff_script)
    fout.close()
    os.system('qsub -P vet_roslin -l h_rt=05:00:00 %s' %(cuff_fname))    

