# -----------------------------------------------------------------------------
# Match porcine Affy probes to KEGG pathways
# -----------------------------------------------------------------------------

#
# This is done by querying Bioconductor porcine.db and org.Ss.eg.db. However,
# it is done thorugh findAttractors in the package attract.

library(attract)
library(RODBC)
library(porcine.db)

## Extract this info from the incidence matrix produced by findAttractors:
probe2kegg<- findAttractors(loring.eset, "time_point", nperm= 1,
    annotation= "porcine.db")@incidenceMatrix
kegg_ids<- rownames(probe2kegg)
probe2kegg_db<- as.data.frame(cbind(
    kegg_id= rep(kegg_ids, times= ncol(probe2kegg)),
    probe_set_id= rep(colnames(probe2kegg), each= length(kegg_ids)),
    present= as.vector(probe2kegg),
    stringsAsFactors= FALSE
))
probe2kegg_db<- probe2kegg_db[which(probe2kegg_db$present == 1), 1:2]
dim(probe2kegg_db)

## Export to postgres
conn<- odbcConnect(dsn= 'pgVitelleschi')
sqlSave(conn, probe2kegg_db, 'probe_to_kegg', rownames= FALSE)
sqlQuery(conn, "comment on table probe_to_kegg is 'Mapping of Affymetrix probes from porcine chip 24K to KEGG pathways. See 20110311_affy_probes_to_kegg.R'")
sqlQuery(conn, 'alter table probe_to_kegg set schema to affymetrix')
odbcClose(conn)