library(biomaRt)
library(RODBC)

mart <- useMart("ensembl")
datasets <- listDatasets(mart)
mart<-useDataset("sscrofa_gene_ensembl", mart)

## For ease of reading: Write the list of attributes to clipboard (paste this in Excel)
writeClipboard(listAttributes(mart)[1][,])
listFilters(mart)

affy<- getBM(attributes= "affy_porcine", mart= mart) ## Values to be returned (all)
affy_annotate<- getBM(
  attributes=c(
    "affy_porcine",
    "ensembl_gene_id",
    "ensembl_transcript_id",
    "external_gene_id"
    ),
  filters= "affy_porcine",
  values= affy[,],
  mart=mart)
dim(affy_annotate)

affy_annotate[affy_annotate == ''] <- NA

colnames(affy_annotate)<- c("probe_set_id", "gene_id", "transcript_id", "gene_symbol")

# -------------------------[ Write to database ]------------------------


conn_db<- odbcConnect(dsn= 'pgVitelleschi', case= 'nochange')
sqlSave(conn_db, affy_annotate, tablename= "porcine_na30_ensembl", rownames= F)
sqlQuery(conn_db, "COMMENT ON TABLE porcine_na30_ensembl IS 'Ensembl annotation of The Affymetrix porcine chip Porcine na30. Done via biomaRt, script 20100711_biomaRt_affy.R'")









