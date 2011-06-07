#______________________________________________________________________________
#
#                        20091021_Biomart_ncbi37_ss9_transcr.R
# 
# Retrieve human pig transcriptome from biomaRt. 
#______________________________________________________________________________

library(RCurl)
library(biomaRt)

#----------------[ Connect to database and retrieve dataset ]------------------

ss9<- useMart(biomart= "ensembl", dataset="sscrofa_gene_ensembl")
h37<- useMart(biomart= "ensembl", dataset="hsapiens_gene_ensembl")

#--------------------------[ Get transcriptomes ]------------------------------

ss9tome<- getBM(attributes= c("ensembl_transcript_id", "chromosome_name", 
  "strand", "transcript_start", "transcript_end"), mart= ss9)

h37tome<- getBM(attributes= c("ensembl_transcript_id", "chromosome_name", 
  "strand", "transcript_start", "transcript_end"), mart= h37)
  
write.table(file="/nfs_netapp/dberaldi/Tritume/ss9transcriptome.txt", ss9tome, 
  quote=FALSE, row.names=FALSE, sep="\t")

write.table(file="/nfs_netapp/dberaldi/Tritume/h37transcriptome.txt", h37tome, 
  quote=FALSE, row.names=FALSE, sep="\t")  

### ----------[ Longer code with - sort of tutorial to biomaRt ]---------------

# listMarts()                 ## See what databases are available
# martdb<- useMart("ensembl") ## Chose a db to work with
# listDatasets(mart=martdb)   ## What datasets are avialable within selected db

### Connect to database and retrieve dataset
# btau<- useMart(biomart= "ensembl", dataset="btaurus_gene_ensembl")

### Show attributes for btau
# write.table(file="/nfs_netapp/dberaldi//Tritume/attributes.txt", 
#    listAttributes(btau), sep="\t", quote=F, row.names=F)

### Get list cow transcriptome with attributes
# cowgenes<- getBM(attributes= c("ensembl_gene_id", "ensembl_transcript_id",
#     "transcript_start", "transcript_end", "start_position", "end_position",
#     "strand"), mart= btau)
