# -----------------------------------------------------------------------------
# Microarray analysis of pig BMDM arrays using LIMMA
# -----------------------------------------------------------------------------

setwd("F:/data/20110204_limma_affymetrix_pig")

library(limma)

# -------------------------[ Import original dataset ]-------------------------

pig_array_ori<- read.table('U:/Documents/affimetrix/lwl_bmdm_lps.txt', 
  sep='\t', header=TRUE, stringsAsFactors= FALSE)

## Input should look like this:
## --------------------------------------------------------------
## Probe_Set_ID	      LWL1a_0	LWL1a_LPS2	LWL1a_LPS7	 ...
## AFFX-BioB-3_at:bioB	304.2370405	362.4969384	361.275347   ...
## AFFX-BioB-5_at:bioB	279.2675236	329.9260767	364.5950085  ...
## AFFX-BioB-M_at:bioB	333.5175803	440.3708	459.5566804  ...
## ...                  ...         ...         ...
## --------------------------------------------------------------

## Remove gene name from probe_set_id to return to the original Affymetrix names
## It is assumed the probe name doesn't have the colon ':' character.
probe_set_id<- sapply(pig_array_ori[,1], function(x) strsplit(x, ':')[[1]][1])

## Note: First take log2 of the intensity and then average
pig_array_log<- data.frame(probe_set_id= probe_set_id, log2(pig_array_ori[, 2:ncol(pig_array_ori)]))

## Average replicates
limma_affy<- cbind(
    LWL1_0h=  rowMeans(pig_array_log[, c(2, 6)]),      
    LWL1_2h=  rowMeans(pig_array_log[, c(3, 7)]),
    LWL1_7h=  rowMeans(pig_array_log[, c(4, 8)]),
    LWL1_24h= rowMeans(pig_array_log[, c(5, 9)]),
    LWL4_0h=  rowMeans(pig_array_log[, c(10, 14)]),
    LWL4_2h=  rowMeans(pig_array_log[, c(11, 15)]),
    LWL4_7h=  rowMeans(pig_array_log[, c(12, 16)]),
    LWL4_24h= rowMeans(pig_array_log[, c(13, 17)]),
    LWL5_0h=  rowMeans(pig_array_log[, c(18, 22)]),
    LWL5_2h=  rowMeans(pig_array_log[, c(19, 23)]),
    LWL5_7h=  rowMeans(pig_array_log[, c(20, 24)]),
    LWL5_24h= pig_array_log[, c(21)]    ## Note: Array in column 25 "s024_LWL5_24h_repl_2" should be excluded (see LabBook 09/07/2010)
    )
rownames(limma_affy)<- probe_set_id
limma_affy[1:10,]


# -----------------------------------------------------------------------------
# Import target file
# -----------------------------------------------------------------------------

# NB: It is essential that the arrays in the *target* dataframe are in the same order as in the inetnsity *matrix* <=== 
targets_ori<- read.table('U:/Documents/affimetrix/pig_bmdm_array_design_avg.txt', sep='\t',header=T)

## Select targets at time points:
for (i in c(2, 7, 24)){
    time_1<- 0
    time_2<- i
    
    tp_compare<- which(targets_ori$time_point == time_1 | targets_ori$time_point == time_2)
    targets<- targets_ori[tp_compare, ]
    
    # Subset to have only required timne points
    limma_affy_sub<- limma_affy[, tp_compare] 
    colnames(limma_affy_sub)<- paste(targets$pig, targets$time_point, sep='_')
    limma_affy_sub[1:10,]
    
    # Specify design matrix (see page 38 'Limma user's guide', para '8.3 Paired Samples')
    pig<- factor(targets$pig)
    time_p<- factor(targets$time_point, levels= unique(targets$time_point))
    
    design<- model.matrix(~ 1 + pig + time_p)
    fit <- lmFit(limma_affy_sub, design)
    fit <- eBayes(fit)
    
    time_p_ind<- grep( 'time_p', names(data.frame(fit$coefficients[1:10,])) ) ## Index of the time point tested
    time_p_name<- names(data.frame(fit$coefficients[1:10,]))[time_p_ind]
    
    topTable(fit, coef= time_p_name)  ##
    de_table<- topTable(fit, coef= time_p_name, number= nrow(fit$coefficients), adjust="BH")
    
    comparison<- paste(time_1, 'vs', time_2) ## Make a nice name for this comparison
    de_table$comparison<- comparison
    
    # Write out to file: Data and graphs
    write.table(de_table, 'lwl_bmdm_lps_limma_toptags.txt', sep='\t', row.names=F, append=T, col.names=F)
    
    tiff(paste(time_p_name, '.tiff', sep=''), res= 150, pointsize= 12, units= 'cm', width = 16, height = 9 )
    par(mfrow=c(1,2), cex=0.8, mar= c(5,4,0,1), oma= c(0,0,4,0))
      hist(de_table$P.Value, main= '', xlab='p-value')
      plot(de_table$AveExpr, de_table$logFC, , col= ifelse(de_table$adj.P.Val<0.01, 'blue', 'black'), pch=16, cex=0.2, 
        main= '', xlab= 'Log average expression', ylab= 'log2 fold change')
      title(outer= T, main= paste("Expression differences at", comparison) )
    dev.off()
}
