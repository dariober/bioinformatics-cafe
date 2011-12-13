## Mapping CAGE library hg95ctrls against S. scrofa genome Ss9.53
## Using bowtie 0.11.3
## -f : Input file is fasta format
## -m 10 : Reject tags which align in more than 10 (-m 10) locations
## -k 10 : For each tag, report up to 10 locations (i.e all of those allowed by -m)

/exports/work/vet_roslin_nextgen/dario/bowtie/bowtie-0.11.3/bowtie \
    Sscrofa9.53.bowtie-index                                       \
    -f /exports/work/vet_roslin_nextgen/dario/fastq/hg95ctrls.fa   \
    -m 10                                                          \
    -k 10                                                          \
    /exports/work/vet_roslin_nextgen/dario/bowtie/output/bwt_hg95ctrls_Ss9_m10.map