#
# Get from biomaRt the annotation for the human refseq_id matching
# the pig affymetrix probes. Matching of Affy probes to human refseq provided
# by Chris Tuggle.
#

library(biomaRt)
library(RODBC)

conn<- odbcConnect(dsn= 'pgVitelleschi')
affy_human<- sqlQuery(conn, 'select distinct refseq_dna_id from affyarray_human_annotation', stringsAsFactors= FALSE)[,1]

affy_human[1:10]
length(affy_human)

mart <- useMart("ensembl")
datasets <- listDatasets(mart)
mart<-useDataset("hsapiens_gene_ensembl", mart)

## For ease of reading: Write the list of attributes to clipboard (paste this in Excel)
writeClipboard(listAttributes(mart)[1][,])
listFilters(mart)

annotate<- getBM(
  attributes=c(
    "refseq_dna",
    "ensembl_gene_id",
    "external_gene_id",
    "description"
    ),
  ## Note: with this filter only NM_ accessions are retrieved
  filters= "refseq_dna",
  values= affy_human,
  mart=mart)
dim(annotate)

sqlSave(conn, annotate)
sqlQuery(conn,"
  update affyarray_human_annotation set
    gene_id = ensembl_gene_id,
    gene_symbol = external_gene_id,
    description = annotate.description
  from annotate where refseq_dna_id = refseq_dna;
  ")
sqlDrop(conn, 'annotate')









