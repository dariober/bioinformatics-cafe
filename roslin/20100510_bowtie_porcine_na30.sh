#
#  Align affymetrix oligonucleotides from chip Procine na30 against
#  S.scrofa transcriptome
#  Report up to 10 valid hits/oligo in order best to worse and grouped by error stratum
#

cd /exports/work/vet_roslin_nextgen/dario/bowtie/output/20100510_porcine_na30

/exports/work/vet_roslin_nextgen/dario/bowtie/current/bowtie \
    /exports/work/vet_roslin_nextgen/dario/bowtie/indexes/Sus_scrofa.Sscrofa9.56.cdna.all \
    -f /exports/work/vet_roslin_nextgen/dario/fasta/affymetrix_porcine_na30/Porcine_na30.fa \
    -a \
    -m 10 \
    -v 2 \
    --best \
    --strata \
    /exports/work/vet_roslin_nextgen/dario/bowtie/output/20100510_porcine_na30/porcine_na30_Sus_scrofa.Sscrofa9.56.cdna.all.bowtiemap
