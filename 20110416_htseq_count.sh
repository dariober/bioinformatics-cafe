#!/bin/bash

#
# Count reads mapping to a GTF feature.
# Input: bam file (typically from tophat)
# Output: Output of HTSeq
#
## qsub -P vet_roslin -l h_rt=06:00:00 job_20110416_htseq_count.sh

source ~/.bash_profile
module add python/2.6.3

gtf='/exports/work/vet_roslin_nextgen/dario/ensembl/release-60/gtf/sus_scrofa/Sus_scrofa.Sscrofa9.60.gtf'

cd /exports/work/vet_roslin_nextgen/dario/tophat/output/20110202_rnaseq_bmdm_lps
python -m HTSeq.scripts.count -q -a 15 -m union -s no -t exon -i gene_id accepted_hits.qnamesorted.sam $gtf > bmdm_lps.htseq

exit
