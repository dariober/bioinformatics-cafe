## Mapping CAGE library hg95ctrls against S. scrofa genome Ss9.53
## Using bowtie 0.11.3
## Keep tags mapping in at most 11 different locations

/exports/work/vet_roslin_nextgen/dario/bowtie/bowtie-0.11.3/bowtie \
    Sscrofa9.53.bowtie-index                                       \
    -f /exports/work/vet_roslin_nextgen/dario/fastq/hg95ctrls.fa   \
    -k 11                                                          \
    --best                                                         \
    /exports/work/vet_roslin_nextgen/dario/bowtie/output/bwt_hg95ctrls_Ss9_k11.map