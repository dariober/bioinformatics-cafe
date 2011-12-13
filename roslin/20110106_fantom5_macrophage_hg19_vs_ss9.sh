#!/bin/bash
PATH=$PATH:/exports/work/vet_roslin_nextgen/dario/bwa/bwa-0.5.8c
cd /exports/work/vet_roslin_nextgen/dario/bwa/output/20110106_fantom5_macro_monocyte_vs_ss9

#bwa aln \
#    -t 4 \
#    /exports/work/vet_roslin_nextgen/dario/bwa/indexes/Sus_scrofa.Sscrofa9.56.dna.toplevel.fa \
#    Macrophage%20-%20monocyte%20derived%20donor1.CNhs10861.11232-116C8.hg19.ctss.fq.gz > 20110106_fantom5_macrophage_hg19_vs_ss9.bwa

bwa samse \
    /exports/work/vet_roslin_nextgen/dario/bwa/indexes/Sus_scrofa.Sscrofa9.56.dna.toplevel.fa \
    20110106_fantom5_macrophage_hg19_vs_ss9.bwa \
    Macrophage%20-%20monocyte%20derived%20donor1.CNhs10861.11232-116C8.hg19.ctss.fq.gz > 20110106_fantom5_macrophage_hg19_vs_ss9.sam


