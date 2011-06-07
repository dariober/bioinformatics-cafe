## CAGE library Test8, filtered to have each read with no base call <20 
## aligned against Sscrofa 9.56 plus human rDNA
## Using bowtie 0.11.3
## Allow up to 3 mismatches in the whole read (SOAP like)
## -q : Input file is fastq format
## -m 10 : Reject tags which align in more than 10 (-m 10) locations
## -k 10 : For each tag, report up to 10 locations (i.e all of those allowed by -m)

/exports/work/vet_roslin_nextgen/dario/bowtie/bowtie-0.11.3/bowtie \
    Sscrofa9.56_plus_Human_ribosomal_DNA_complete_repeating_unit                                       \
    -q /exports/work/vet_roslin_nextgen/dario/fastq/test8/s_6_test8_qs20.fastq \
    -m 10                                                          \
    -k 10                                                          \
    -v 3                                                           \
    /exports/work/vet_roslin_nextgen/dario/bowtie/output/bwt_Test8_ss9_plus_human_rDNA.map