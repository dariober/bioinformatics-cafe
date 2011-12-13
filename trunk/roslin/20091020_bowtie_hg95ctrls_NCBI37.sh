## Mapping CAGE library hg95ctrls against H. sapiens genome NCBI37 (hg19)
## Using bowtie 0.11.3

/exports/work/vet_roslin_nextgen/dario/bowtie/bowtie-0.11.3/bowtie \
    h_sapiens_37_asm                                               \
    -f /exports/work/vet_roslin_nextgen/dario/fastq/hg95ctrls.fa   \
    -k 1                                                           \
    --best                                                         \
    /exports/work/vet_roslin_nextgen/dario/bowtie/output/bwt_hg95ctrls_NCBI37.map