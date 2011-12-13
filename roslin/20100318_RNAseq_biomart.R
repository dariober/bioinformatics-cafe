library(RODBC)

conn_db<- odbcConnect(dsn= "pgCAGE", case= "nochange", uid= "dberaldi", pwd= "mypassword")
fpkm<- sqlQuery(conn_db,
  "select fpkm_ctrl.*, fpkm_lps.fpkm_lps,  log(2, (fpkm_lps.fpkm_lps/fpkm_ctrl.fpkm_ctrl)::numeric) AS fold_change from 
    (select distinct transcript_id, fpkm AS fpkm_ctrl from cufflinks_transcript_gtf where source like '20100317_RNAseq_CTRL') AS fpkm_ctrl 
    INNER JOIN 
    (select distinct transcript_id, fpkm AS fpkm_lps from cufflinks_transcript_gtf where source like '20100317_RNAseq_LPS') AS fpkm_lps ON
    fpkm_ctrl.transcript_id= fpkm_lps.transcript_id
    order by fold_change;"
  )
dim(fpkm)
fpkm[1:10,]

tnf<- fpkm[fpkm$transcript_id=='ENSSSCT00000001533',]
hprt<- fpkm[fpkm$transcript_id=='ENSSSCT00000013865',]
actb<- fpkm[fpkm$transcript_id=='ENSSSCT00000008324',]

par(mfrow=c(1,2), cex= 0.7)
plot(fpkm$fold_change, ylab= "Fold change - log2(lps/ctrl)", main= "CTRL vs LPS\nFPKM differences", xlab="Transcript ID")

## abline(h= 0, col= "black", lty= "dotted")
abline(h=tnf$fold_change, lty="dotted", col="red")
abline(h=hprt$fold_change, lty="dotted", col="blue")
abline(h=actb$fold_change, lty="dotted", col="blue")
text(labels= "TNF-alpha", y=tnf$fold_change+0.5, x=0, adj=0)
text(labels= "HPRT", y=hprt$fold_change-0.5, x= 8000, adj=0)
text(labels= "ACTB", y=actb$fold_change+0.5, x= 0, adj=0)

##hist(fpkm$fold_change, freq= FALSE, ylim= c(0,0.8))
d= density(fpkm$fold_change)
plot(d)

plot(pnorm(fpkm$fold_change))

dim(fpkm[fpkm$fold_change <= -tnf$fold_change, ])

hist(fpkm$fold_change)

zscore<- function(x){
    meanx<- mean(x)
    sdx<- sd(x)
    return((x - meanx) / sdx)
    }
z_fold<- zscore(fpkm$fold_change)
pval<- 2 * (1 - abs(pnorm(q= z_fold, mean= mean(z_fold), sd= sd(z_fold))))
plot(z_fold)

odbcClose(conn_db)

# ----------------------[ Tritume ]-------------------
## BiomaRt doesn't work!
sscrofa<- useMart(biomart= "ensembl", dataset="sscrofa_gene_ensembl")
listAttributes(sscrofa)
rnaseq_ann<- getBM(mart= sscrofa, filters='ensembl_transcript_id', 
  values= 'ENSSSCT00000001533',
  attributes = c("entrezgene")
  )
