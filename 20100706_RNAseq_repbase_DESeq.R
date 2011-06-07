library(DESeq)

## This table produced by 20100706_sam_bowtie_rebase.sql
## See also LabBook 06/07/2010

countsTable<- read.table('D:/Tritume/RNAseq_repbase.txt', sep= '\t', header= T, comment.char= '', stringsAsFactors= FALSE)

## Exclude simple repeats
no_simple<- (1:nrow(countsTable) %in% grep('.*Simple_repeat.*', countsTable$rname)) == FALSE

## Exclude tRNA
no_tRNA<- (1:nrow(countsTable) %in% grep('^tRNA', countsTable$rname)) == FALSE

## New dataset without some families of repeats
countsTable<- countsTable[no_simple & no_tRNA, ]
dim(countsTable)

names(countsTable)[1]<- 'gene' ## Just to be consistent with the DESeq vignettes
countsTable[is.na(countsTable)]<- 0

rownames(countsTable)<- countsTable$gene
countsTable<- countsTable[,-1]

countsTable[1:10,]

conds<- c('CTRL', 'LPS')

cds<- ?newCountDataSet(countsTable, conds)
head(counts(cds))

## The size factor is value associated to each librairy to scale all the libraries to the same size
cds<- estimateSizeFactors(cds)
sizeFactors(cds)

## Variance estimation
cds<- estimateVarianceFunctions(cds, pool= TRUE)

# ----------------------[ Some diagnostics... Possibly Not Applicable here ]------------------
scvPlot(cds)
diagForT<- varianceFitDiagnostics(cds, 'CTRL')
head(diagForT)

smoothScatter( log10(diagForT$baseMean), log10(diagForT$baseVar) )
lines( log10(fittedBaseVar) ~ log10(baseMean),
 diagForT[ order(diagForT$baseMean), ], col="red" )

# -------------------------[ Estimate differential expression ]-------------------------------

res<- nbinomTest(cds, 'LPS', 'CTRL')

plot(res$baseMean, res$log2FoldChange, log="x", pch=20, cex= 1,plot
 col = ifelse( res$padj < .1, "red", "black" ), xlab= 'Base mean expression', ylab= 'Log2 fold change')
savePlot('M:/Documents/LabBook/LabBook_Figures/20100706_DESeq_repbase.emf', 'emf')

res[1:10,]
res[res$padj < 0.1, ]

