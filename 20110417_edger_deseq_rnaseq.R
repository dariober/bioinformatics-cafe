# ------------------------------------------------------------------------------
# edgeR differential expression RNAseq
# ------------------------------------------------------------------------------

library(edgeR)
library(RODBC)
library(DESeq)

conn<- odbcConnect(dsn= 'pgVitelleschi')
htseq<- sqlQuery(conn, 'SELECT * FROM htseq_count')
htseq<- reshape(htseq, v.names= 'count', idvar= 'gene_id', timevar= 'source', direction= 'wide')
colnames(htseq)<- c('gene_id', 'am_ctrl', 'am_lps', 'bmdm_ctrl', 'bmdm_lps')
rownames(htseq)<- htseq$gene_id
htseq<- htseq[, -1]
htseq[is.na(htseq)]<- 0
htseq[1:10,]
#                    am_ctrl am_lps bmdm_ctrl bmdm_lps
# ENSSSCG00000000002       5      4       324      117
# ENSSSCG00000000003     296    187       203      163
# ENSSSCG00000000004       5      2         5        0
# ENSSSCG00000000005      54     70       103       58
# ENSSSCG00000000006       3      1        47        6


# ------------------------------------------------------------------------------
# DESeq: AM CTRL vs AM LPS and BMDM ctrl vs lps
# ------------------------------------------------------------------------------

## AM
raw.am<- htseq[, c(1,2)]
group<- c('ctrl', 'lps')

cds <- newCountDataSet( raw.am, group )
cds <- estimateSizeFactors( cds )       ## Factor to adjust library sizes to account for differences in seq depth
## method= 'blind' allows for no replication.
cds <- estimateVarianceFunctions( cds, method= 'blind' )
res.am <- nbinomTest( cds, "ctrl", "lps")  ## Binomial test for D.E.
head(res)

## BMDM
raw.bmdm<- htseq[, c(3,4)]
group<- c('ctrl', 'lps')

cds <- newCountDataSet( raw.bmdm, group )
cds <- estimateSizeFactors( cds )       ## Factor to adjust library sizes to account for differences in seq depth
## method= 'blind' allows for no replication.
cds <- estimateVarianceFunctions( cds, method= 'blind' )
res.bmdm<- nbinomTest( cds, "ctrl", "lps")  ## Binomial test for D.E.

## Plots:
tiff('M:/Documents/LabBook/LabBook_Figures/20110417_deseq_diffexpr.tiff', width= 16, height= 7, res= 150, units= 'cm', pointsize= 9)
par(mfrow= c(1,2), cex= 0.85, las= 1, mgp= c(2,1,0), mar= c(3.5, 3, 2, 0.1), oma= c(0,0,0,2), xpd= TRUE)
  with(res.am, {
    plot(baseMean, log2FoldChange, log= 'x', pch=20, cex=.1, col = ifelse( padj < .1, "red", "black" ), main= 'Alveol. Macroph.', xlab= 'Base mean', ylab= 'Log2 fold change')
  })
  abline(h= c(-2,2), col= 'dodgerblue', lwd= 1.25, xpd= FALSE)
  with(res.bmdm, {
    plot(baseMean, log2FoldChange, xaxt= 'n', log= 'x', pch=20, cex=.1, col = ifelse( padj < .1, "red", "black" ), main= 'Bone Marrow Der. Macroph.', xlab= 'Base mean', ylab= '')
    x<- axis(side= 1, labels= FALSE)
    axis(side= 1, at= x, labels= formatC(x, format= 'd'))
  })
  abline(h= c(-2,2), col= 'dodgerblue', lwd= 1.25, xpd= FALSE)
dev.off()
graphics.off()

## Export:
res.am$comparison<- 'am'
res.bmdm$comparison<- 'bmdm'
head(res.am)
res<- rbind(res.am, res.bmdm)
# sqlQuery(conn, 'DROP TABLE deseq_nbinom')
sqlSave(conn, res, 'deseq_nbinom', rownames= FALSE)
sqlQuery(conn, 'ALTER TABLE deseq_nbinom SET SCHEMA "Pigs"')
sqlQuery(conn, "COMMENT ON TABLE deseq_nbinom IS 'Output of DESeq function nbinomTest on AM and BMDM rnaseq. Alignment from 20110202. See 20110417_edger_deseq_rnaseq.R. Labbook 16/04/2011'")

# ------------------------------------------------------------------------------
# NOT USED:
# EdgeR: AM CTRL vs AM LPS
# ------------------------------------------------------------------------------

raw.data<- htseq[, c(3,4)]
group<- c('ctrl', 'lps')

d <- DGEList(counts = raw.data, group = group)
d <- estimateCommonDisp(d)

de.com <- exactTest(d)
names(de.com)
names(de.com$table)
topTags(de.com)

## See raw counts for the top d.e. genes
detags.com <- rownames(topTags(de.com)$table) 
d$counts[detags.com, ]
 
## Count features with p<0.01
n.de<- sum(de.com$table$p.value < 0.01) 
sum(p.adjust(de.com$table$p.value, method = "BH") < 0.01) ## Adjust for multiple testing

de.table<- topTags(de.com, n= nrow(de.com$table))$table
length(de.table$FDR[de.table$FDR < 0.01])

windows(width= 10/2.54, height= 10/2.54)
par(cex= 0.85, las= 1)
plot(de.table$logConc, de.table$logFC, pch= 20, col= ifelse(de.table$FDR < 0.01, 'red', 'grey40'))
graphics.off()

# -----------------------------------------------------------------------------
# Comparison in log2fc between RNAseq and Affymetrix arrays
# -----------------------------------------------------------------------------

library("geneplotter")  ## from BioConductor
require("RColorBrewer") ## from CRAN

## This table produced by 20110418_affy_rnaseq_bmdm.sql.
## Affymetrix probes annotated to gene_id through C.Tuggle annotation (gene_name)
## gene_name linked to ensembl_gene_id via biomart. Only expression for LWL1
rnaffy<- read.table('D:/Tritume/affy_rnaseq.txt', sep= '\t', header=TRUE)
dim(rnaffy) ## 7068

## Remove genes with Infintiy log2fc in rnaseq
rnaffy<- rnaffy[abs(rnaffy$log2fc_rnaseq) < 20,]
dim(rnaffy) ## 6668

x<- rnaffy$log2fc_affy
y<- rnaffy$log2fc_rnaseq
mylm<- lm(y ~ x + I(x^2) + I(x^3) + I(x^4)+ I(x^5))
summary(mylm)
anova(mylm)
newx<- seq(-5, 8) 
prd<- predict(mylm, newdata= data.frame(x=newx), interval = c("prediction"), type="response")

tiff('M:/Documents/LabBook/LabBook_Figures/20110417_affy_vs_rnaseq.tif', res= 150, pointsize= 8, units= 'cm',
    width = 7, height = 7)
par(las= 1, cex= 0.85, mgp=c(2.25, 0.75, 0))
smoothScatter(rnaffy$log2fc_affy, rnaffy$log2fc_rnaseq, bandwidth=0.001, nrpoints= 200, main= "Affymetrix vs RNAseq\n BMDM - LWL1", yaxt='n',
              xlab= 'Affymetrix log2 fc', ylab= 'RNAseq log2 fc')
yx<- axis(side= 2, at= seq(-6, 8, by= 2))
lines(newx, prd[,1], col="red", lty=2, lwd= 1)
lines(newx, prd[,2], lty=2, lwd= 1, col= 'grey')
lines(newx, prd[,3], lty=2, lwd= 1, col= 'grey')
abline(h= yx, v= seq(-2,6, by=2), col= 'grey60', lty= 'dotted')
dev.off()

cor.test(x,y)
summary(lm1)


#
# TRITUME
#

plot(x,y)
scatter.smooth(x, y, col= 'blue', lwd= 2)