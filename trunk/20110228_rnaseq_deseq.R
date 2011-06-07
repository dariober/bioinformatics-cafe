#
#  RNAseq analysis of differential expression using DESeq
#

library(RODBC)
library(DESeq)
library(limma)
library(gplots)

# ------------------------------------------------------------------------------
# Import and prepare data
# ------------------------------------------------------------------------------

conn<- odbcConnect(dsn= 'pgVitelleschi')

sqlQuery(conn,
    "select cross_tab($$ select feature, source, count from htseq_count_gene_id $$, 'htseq_count_gene_id_ct');"
    )
htseq_counts<- sqlFetch(conn, 'htseq_count_gene_id_ct')
sqlQuery(conn, "drop table htseq_count_gene_id_ct");
odbcClose(conn)

head(htseq_counts)

htseq_counts[is.na(htseq_counts)] <- 0
gene_id<- htseq_counts[,1]

countsTable<- htseq_counts[,-1]
rownames(countsTable)<- gene_id
countsTable<- countsTable[rowSums(countsTable)>0, ]
countsTable[1:100,]
dim(countsTable)

# V
vcounts<- countsTable
vcounts[vcounts > 0]<- TRUE

rnaseq_venn<- venn(vcounts)

rnaseq_venn<- vennCounts(countsTable)
vennDiagram(rnaseq_venn)

# ------------------------------------------------------------------------------
# Conditions: CTRL / LPS
# ------------------------------------------------------------------------------

## Make sure columns in expression matrix match conditions vector:
## conds<- as.factor(c('CTRL', 'LPS', 'CTRL', 'LPS'))
conds<- as.factor(c('CTRL', 'LPS'))

# ------------------------------------------------------------------------------
# Differential expression
# ------------------------------------------------------------------------------

cds <- newCountDataSet( countsTable, conds )
cds <- estimateSizeFactors( cds )       ## Factor to adjust library sizes to account for differences in seq depth
cds <- estimateVarianceFunctions( cds, method= 'blind' ) ## For each gene, use the mean across all samples to estimate the associated variance

res <- nbinomTest( cds, "CTRL", "LPS")  ## Binomial test for D.E.
head(res)

plot( 
   res$baseMean,
   res$log2FoldChange, 
   log="x", pch=20, cex=.1, 
   col = ifelse( res$padj < .1, "red", "black" ) )

resSig <- res[ res$padj < 0.1, ]
dim(resSig)

# ------------------------------------------------------------------------------
# Diagnostics
# ------------------------------------------------------------------------------

windows(width= 10/2.54, height= 10/2.54)
par( fg= 'grey20' )
par(
    mar= c(4, 4, 2, 1),
    pch= 19,               ## solid circle point
    col= 'grey20',
    col.lab= 'grey40',
    col.axis= 'grey40', 
    las= 1, 
    cex= 0.95,
    mgp= c(2.2, 0.75, 0), ## The margin line (in mex units) for the axis title, axis labels and axis line
    tcl = -0.25,           ## Length of tickmarks
    bty='l')               ## Plot only x and y axis

scvPlot(cds, ylim= c(0, 3))

## Does the base variance follow the empirical variance?
diagForT<- varianceFitDiagnostics( cds, "CTRL" )
head( diagForT )

smoothScatter( log10(diagForT$baseMean), log10(diagForT$baseVar),
    xlab= 'Base mean', ylab= 'Base variance', main= 'Variance vs mean' )
lines( log10(fittedBaseVar) ~ log10(baseMean),
    diagForT[ order(diagForT$baseMean), ], col="red" )

par( mfrow=c(1,2 ) )
residualsEcdfPlot( cds, "CTRL" )
residualsEcdfPlot( cds, "LPS" )
hist(diagForT$pchisq)

# ------------------------------------------------------------------------------
# Comparison RNAseq - affymetrix
# ------------------------------------------------------------------------------

conn<- odbcConnect(dsn= 'pgVitelleschi')

# --------------------------------[ TRITUME ]-----------------------------------

hsb2<-read.table("http://www.ats.ucla.edu/stat/R/notes/hsb2.csv", sep=',', header=T)
attach(hsb2)

hw<-(write>=60)
hm<-(math >=60)
hr<-(read >=60)
c3<-cbind(hw, hm, hr)

vennCounts
