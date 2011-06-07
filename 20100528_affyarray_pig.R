#
# Microarray
#

library(limma)
# ----------------------------------[ Limma ]----------------------------------

pig_array_ori<- read.table('U:/Documents/affimetrix/lwl_bmdm_lps.txt', sep='\t', header=TRUE, stringsAsFactors= FALSE)

## Input should look like this:
## --------------------------------------------------------------
## Probe_Set_ID	      LWL1a_0	LWL1a_LPS2	LWL1a_LPS7	 ...
## AFFX-BioB-3_at:bioB	304.2370405	362.4969384	361.275347   ...
## AFFX-BioB-5_at:bioB	279.2675236	329.9260767	364.5950085  ...
## AFFX-BioB-M_at:bioB	333.5175803	440.3708	459.5566804  ...
## ...                  ...         ...         ...
## --------------------------------------------------------------

colnames(pig_array_ori)

## Remove gene name from probe_set_id to return to the original Affymetrix names
## It is assumed the probe name doesn't have the colon ':' character.
probe_set_id<- sapply(pig_array_ori[,1], function(x) strsplit(x, ':')[[1]][1])

## Note: First take log2 of the intensity and then average
limma_affy<- data.frame(probe_set_id= probe_set_id, log2(pig_array_ori[, 2:ncol(pig_array_ori)]))

## Average replicates
limma_affy<- cbind(limma_affy[,1],              
    rowMeans(limma_affy[, c(2, 6)]),      
    rowMeans(limma_affy[, c(3, 7)]),
    rowMeans(limma_affy[, c(4, 8)]),
    rowMeans(limma_affy[, c(5, 9)]),
    rowMeans(limma_affy[, c(10, 14)]),
    rowMeans(limma_affy[, c(11, 15)]),
    rowMeans(limma_affy[, c(12, 16)]),
    rowMeans(limma_affy[, c(13, 17)]),
    rowMeans(limma_affy[, c(18, 22)]),
    rowMeans(limma_affy[, c(19, 23)]),
    rowMeans(limma_affy[, c(20, 24)]),
    rowMeans(limma_affy[, c(21, 25)])
    )
## Remove column of probe names, set rownames to probe names
rownames(limma_affy)<- probe_set_id
limma_affy<- limma_affy[, 2:ncol(limma_affy)]

## ===> NB: It is essential that the arrays in the *target* dataframe are in the same order as in the inetnsity *matrix* <=== 
##
targets<- read.table('U:/Documents/affimetrix/pig_bmdm_array_design_avg.txt', sep='\t',header=T)

## Select targets at 0h and 7h
h0_7<- which(targets$time_point == 0 | targets$time_point == 7)
targets<- targets[h0_7, ]

# Subset to have only 0 and 7h
limma_affy<- limma_affy[, h0_7] ## Arrays at 0h and 7h
colnames(limma_affy)<- paste(targets$pig, targets$time_point, sep='_')
limma_affy[1:10,]

# Specify design matrix (see page 38 'Limma user's guide', para '8.3 Paired Samples')
pig<- factor(targets$pig)
time_p<- factor(targets$time_point, levels= unique(targets$time_point))

design<- model.matrix(~ 1 + pig + time_p)
fit <- lmFit(limma_affy, design)
fit <- eBayes(fit)
topTable(fit, coef="time_p7")

de_table<- topTable(fit, coef="time_p7", number= nrow(fit$coefficients), adjust="BH")

par(mfrow=c(1,2), cex=0.8)
  hist(de_table$P.Value, main= 'Histogram of p-values for all probes', xlab='p-value')
  plot(de_table$AveExpr, de_table$logFC, , col= ifelse(de_table$adj.P.Val<0.01, 'blue', 'black'), pch=16, cex=0.2, 
    main= 'Differential expression', xlab= 'Log average expression', ylab= 'log2 fold change')

length(de_table$adj.P.Val[de_table$adj.P.Val<0.01 & de_table$logFC > 0])
length(de_table$adj.P.Val[de_table$adj.P.Val<0.01 & de_table$logFC < 0])
  
# Write out to file
write.table(de_table, 'U:/Documents/affimetrix/lwl_bmdm_lps_avg_ct_0h_7h.limma', sep='\t', row.names=F)

# -----------------------------------------------------------------------------
#                                Affymetrix vs RNAseq 
# -----------------------------------------------------------------------------

## Table produced by 20100413_affyarray_RNAseq.sql
affy_rnaseq<- read.table('C:/Tritume/affy_rnaseq_db.txt', sep='\t', header=T)
dim(affy_rnaseq)
targets<- read.table('U:/Documents/affimetrix/pig_bmdm_array_design.txt', sep='\t',header=T)

# ----------------[ RNAseq vs Arrays for some genes of interest ]--------------

time_points<- unique(targets$time_point)
pig_colours<- c('red', 'blue', 'green')
genes<-c('TNF', 'CD40', 'IL1A', 'CSF3', 'CCL4', 'CCL3L1')

par(mfrow=c(3,2), mar=c(1,2.5,1,1), oma=c(2,1.5,2,0), cex= 0.9)
  for (gene in genes){
    max_y<- max(log2(mean(affy_rnaseq[affy_rnaseq$gene_symbol == gene, 10:33]))) ## Make sure your are picking the right cols!
    min_y<- min(log2(mean(affy_rnaseq[affy_rnaseq$gene_symbol == gene, 10:33])))

    fpkm_ctrl<- mean(affy_rnaseq$fpkm_ctrl[affy_rnaseq$gene_symbol == gene], na.rm=T)
    fpkm_ctrl<- ifelse(is.nan(fpkm_ctrl), 0.001, fpkm_ctrl)
    fpkm_lps<- mean(affy_rnaseq$fpkm_lps[affy_rnaseq$gene_symbol == gene], na.rm=T)
    fpkm_lps<- ifelse(is.nan(fpkm_lps), 0.001, fpkm_lps)

    ymax<- max(c(log(fpkm_ctrl), log2(fpkm_lps), max_y))
    barplot(log2(c(fpkm_ctrl, fpkm_lps)), space= c(6,0), axes=F, width= 1, col= 'lightblue', 
        xlim= c(min(time_points), max(time_points)), ylim= c(0, ymax)
        )
    axis(side=2, las=1)
    axis(side=1, at= c(0,2,7,24), labels=F)
    legend('bottomright', bty='n', legend=gene, cex=1.2)
    
    for(i in 1:length(unique(targets$pig))){
        pig<- unique(targets$pig)[i]
        points(time_points, log2(mean(affy_rnaseq[affy_rnaseq$gene_symbol == gene, (which(targets$cell_culture_replicate == 1 & targets$pig == pig) + 9)])), type= 'o', pch=16, col=pig_colours[i])
        points(time_points, log2(mean(affy_rnaseq[affy_rnaseq$gene_symbol == gene, (which(targets$cell_culture_replicate == 2 & targets$pig == pig) + 9)])), type= 'o', pch=16, col=pig_colours[i])
        }   
    }
    title(main= 'Array and RNAseq expression for some genes of interest', outer=TRUE)
    mtext(text='Log2 expression', side= 2, line= 0, cex= 0.9, outer=T)
    mtext(side=1, line= 0, text= rep('Time point (h)', 2) , cex=0.9, outer=TRUE, at=c(0.25, 0.75))
savePlot('M:/Documents/LabBook/LabBook_Figures/20100413_array_vs_rnaseq.emf', 'emf')

# --------------------------[ Overview for all DE genes ]----------------

de_genes<- with(affy_rnaseq,
  # Restrict comparison to transcripts found in both RNAseq CTRL and LPS, FDR(RNAseq)<0.01
  affy_rnaseq[is.na(fpkm_ctrl) == FALSE & is.na(fpkm_lps) == FALSE & is.na(affy_log2fc) == FALSE & rnaseq_fdr < 0.01,]
  )
dim(de_genes)
names(de_genes)

de_avg<- aggregate(data.frame(affy_log2fc= de_genes$affy_log2fc, rnaseq_log2fc= de_genes$rnaseq_log2fc, affy_fdr= de_genes$affy_fdr, rnaseq_fdr= de_genes$rnaseq_fdr), 
    by= list(transcript_id= de_genes$transcript_id), mean)
de_avg<- de_avg[order(de_avg$rnaseq_log2fc), ]

sel_t<- unique(de_genes$transcript_id[de_genes$gene_symbol %in% genes]) ## Transcripts matching genes of interest from above plot

length(de_avg$transcript_id[de_avg$rnaseq_fdr<0.01])
length(de_avg$transcript_id[de_avg$rnaseq_fdr<0.01]) / nrow(de_avg)
length(de_avg$transcript_id[de_avg$affy_fdr<0.05])
length(de_avg$transcript_id[de_avg$affy_fdr<0.05]) / nrow(de_avg)

lm_de1<- lm(de_avg$affy_log2fc ~ de_avg$rnaseq_log2fc)
lm_de2<- lm(de_avg$affy_log2fc ~ de_avg$rnaseq_log2fc + I(de_avg$rnaseq_log2fc^2))

## Plot comparison
## windows(width= 17/2.54, height= 17/2.54)
plot(de_avg$rnaseq_log2fc, de_avg$affy_log2fc, 
  ylim= c(min(de_avg$rnaseq_log2fc, de_avg$affy_log2fc), max(de_avg$rnaseq_log2fc, de_avg$affy_log2fc)),
  xlim= c(min(de_avg$rnaseq_log2fc, de_avg$affy_log2fc), max(de_avg$rnaseq_log2fc, de_avg$affy_log2fc)),
  xlab= 'Log2 RNAseq', ylab= 'Log2 Array', main= 'Fold change for transcripts in RNAseq and Affymetrix', cex.main= 1.25,
  cex= 0.95)
points(de_avg$rnaseq_log2fc[de_avg$affy_fdr>0.05], de_avg$affy_log2fc[de_avg$affy_fdr>0.05], col='blue')
points(de_avg$rnaseq_log2fc, lm_de2$fitted.values, type='l', col= 'red', lwd= 2, lty= 'dotted')
##
## This points() highlights the genes from the previuos plot
##
## points(de_avg$rnaseq_log2fc[de_avg$transcript_id %in% sel_t], 
##     de_avg$affy_log2fc[de_avg$transcript_id %in% sel_t], col= 'red', pch= 16, cex= 0.75)
abline(a= 0, b= 1, lty= 'dotted')
legend('topleft', bty='n', cex= 1.1, legend= paste('Affy probes with p > 0.05 (n= ', length(de_avg$affy_log2fc[de_avg$affy_fdr>=0.05]), ')', sep=''),
    pch= 16, col= 'blue')
legend('bottomright', bty='n', cex= 1.1, legend= paste('n transcripts=', nrow(de_avg)))
savePlot('M:/Documents/LabBook/LabBook_Figures/20100413_affyarray_RNAseq_Fig_1.emf', 'emf')

summary(lm_de2)
summary(lm_de1)
anova(lm_de1, lm_de2)


# -----------------------------------------------------------------------
#                                Tritume
# -----------------------------------------------------------------------

de_genes[1:10,]
affy_rnaseq<- affy_rnaseq[order(affy_rnaseq$rnaseq_log2fc), ]
plot(affy_rnaseq$rnaseq_log2fc, affy_rnaseq$affy_log2fc, xlim=c(1, 6), ylim= c(1,6))
cor.test(affy_rnaseq$rnaseq_log2fc, affy_rnaseq$affy_log2fc)


round(fpkm_labels/max(fpkm_labels),1)

# -----------------------------[ Affymetrix vs RNAseq ]------------------------

affy<- read.table('C:/Tritume/affy_rnaseq.txt', sep='\t', header=TRUE)
dim(affy)

boxplot(affy$cv_intensity) ## Coefficient of variation (with two replicate not very meaningful though!)
affy[1:10,]

# -------------------------[ Home made diff. expr. test ]----------------------

probe<- rownames(limma_affy)
affy_pvalue<- vector(length=length(probe))

time_contr<- which(targets$time_point == 0 | targets$time_point == 7)
for(i in 1:length(probe)){
    avg_intensity<- as.numeric(limma_affy[i,time_contr])
    trans_lm<- lm(avg_intensity ~ as.factor(targets$time_point[time_contr]) + targets$pig[time_contr])
    summ_lm<- summary(trans_lm)
    affy_pvalue[i]<- summ_lm$coefficients[2, 4]
    }

affy_pvalue_2<- vector(length=length(probe))
time_contr<- which(targets$time_point == 0 | targets$time_point == 7)
for(i in 1:length(probe)){
    avg_intensity<- as.numeric(limma_affy[i,time_contr])
    trans_lm<- lm(avg_intensity ~ 0 + as.factor(targets$time_point[time_contr]) + targets$pig[time_contr])
    summ_lm<- summary(trans_lm)
    affy_pvalue_2[i]<- summ_lm$coefficients[2, 4]
    }

length(affy_pvalue[affy_pvalue<0.01])
affy_adjp<- p.adjust(affy_pvalue, method= 'BH')
length(affy_adjp[affy_adjp<0.1])

hist(affy_adjp)
hist(affy_pvalue)

## pvalue<- aggregate(affy$PValue, by= list(affy$transcript_id), min)
## names(pvalue)<- c('transcript_id', 'rnaseq_pvalue')
## pvalue<- data.frame(cbind(pvalue, affy_pvalue= affy_pvalue))

hist(pvalue$affy_pvalue, breaks=40)
length(pvalue$affy_pvalue[pvalue$affy_pvalue<0.01])
plot(pvalue$rnaseq_pvalue, pvalue$affy_pvalue, xlim=c(0, 0.01), ylim=c(0, 0.01))

# --------------------------[ Affy log-fold change ]-----------------------------

lfc_affy<- aggregate(affy$avg_intensity, by= list(affy$transcript_id, affy$time_point), mean)
lfc_affy<- lfc_affy[order(lfc_affy$Group.2, lfc_affy$Group.1), ]
names(lfc_affy)<- c("transcript_id", "time_point", "avg_expr")

lfc_rnaseq<- aggregate(affy$logFC, by= list(affy$transcript_id), min)
lfc_rnaseq<- lfc_rnaseq[order(lfc_rnaseq$Group.1), ]

lfc_affy<- data.frame(transcript_id= unique(lfc_affy$transcript_id), 
  affy_fc= log2(lfc_affy$avg_expr[lfc_affy$time_point == 7] / lfc_affy$avg_expr[lfc_affy$time_point == 0]),
  affy_pvalue= pvalue$affy_pvalue, rnaseq_fc= lfc_rnaseq$x, rnaseq_pvalue= pvalue$rnaseq_pvalue)
lfc_affy[1:10,]

plot(lfc_affy$affy_fc, lfc_affy$rnaseq_fc, xlim=c(-4, 6), ylim=c(-4, 6))



hist(log(affy$avg_intensity))

summ_lm$coefficients[2, 4]

# This would be more appropriate but mpore difficult to interpret. affy from postgresql has cell_culture_replicate NOT grouped
trans_lm<- lme(intensity ~ time_point + pig, random= ~cell_culture_replicate|pig, data= affy[affy$transcript_id == 'ENSSSCT00000000241',], method='ML')

dim(pig_array_n)
pig_array_n[1:10,]

intensity[1:10]
probe_set_id[1:10]
slide_id[1:10]


## File with tissue replicates averaged
pig_array<- read.table('U:/Documents/affimetrix/lwl_bmdm_lps_avg_ct.txt', sep='\t', header=TRUE)
targets<- read.table('U:/Documents/affimetrix/pig_bmdm_array_design_avg.txt', sep='\t',header=T)

# Remove column of probe names and set rownames
limma_affy<- log(pig_array[ ,2:ncol(pig_array)])
rownames(limma_affy)<- pig_array[,1]
limma_affy[1:10,]

# Prepare design matrix
targets$Target<- paste('tp_', targets$time_point, 'h', sep='')
lev<- lev<- unique(targets$Target)
f <- factor(targets$Target, levels=lev)
design <- model.matrix(~0+f)
colnames(design) <- lev

cont_0h_7h<- ?makeContrasts("tp_7h-tp_0h", levels=design)

# Fit linear model
fit<- lmFit(limma_affy, design)
fit <- eBayes(fit)

# Estimate expression differences for contrasts of interest
fit2 <- contrasts.fit(fit, cont_0h_7h)
fit2 <- eBayes(fit2)
# Summarize 
topTable(fit2, number= 10, adjust="BH")
de_table<- topTable(fit2, number= nrow(fit$coefficients), adjust="BH")
length(de_table$adj.P.Val[de_table$adj.P.Val<0.01 & de_table$logFC > 0])
length(de_table$adj.P.Val[de_table$adj.P.Val<0.01 & de_table$logFC < 0])

par(mfrow=c(1,2), cex=0.8)
  hist(de_table$P.Value, main= 'Histogram of p-values for all probes', xlab='p-value')
  plot(de_table$AveExpr, de_table$logFC, , col= ifelse(de_table$adj.P.Val<0.01, 'blue', 'black'), pch=16, cex=0.2, 
    main= 'Differential expression', xlab= 'Log average expression', ylab= 'log2 fold change')
#  plot(de_table$logFC, de_table$B, , col= ifelse(de_table$adj.P.Val<0.01, 'blue', 'black'), pch=16, cex=0.2)
