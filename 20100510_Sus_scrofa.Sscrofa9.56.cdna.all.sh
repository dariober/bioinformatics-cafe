# 
#  Prepare bowtie index for S.scrofa transcriptome 9.56
#  FASTA file downloaded from  
#  wget ftp://ftp.ensembl.org/pub/release-56/fasta/sus_scrofa/cdna/Sus_scrofa.Sscrofa9.56.cdna.all.fa.gz
#

cd /exports/work/vet_roslin_nextgen/dario/bowtie/indexes
/exports/work/vet_roslin_nextgen/dario/bowtie/current/bowtie-build /exports/work/vet_roslin_nextgen/dario/cdna/Sus_scrofa.Sscrofa9.56.cdna.all.fa Sus_scrofa.Sscrofa9.56.cdna.all